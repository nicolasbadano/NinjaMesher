// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#include "Layers.hpp"
#include "Cutter.hpp"

#include <cstdio>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <array>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>


namespace ninja {

namespace {

// --- internal clamp/march constants (defaults, NOT dict knobs --
// kept internal until the benchmark suite shows a case that needs a
// knob) ------------------------------------------------
constexpr int kLaplacianPasses = 3;
// MEASURED (offset_sphere): 0.5 left enough per-point direction
// divergence between neighbours that many side/bottom quads came out
// too twisted for checkMesh's per-face pyramid check (see prismValid's
// twist guard); more neighbour weight smooths that divergence away
// before it ever reaches the clamps.
constexpr double kLaplacianBlend = 0.25;   // new = blend*old + (1-blend)*neighbourAvg
constexpr double kMaxAngleDeg = 75.0;      // max angle between march dir and -frontNormal
// Per-edge tangential stretch/compression bound. MEASURED (offset_sphere):
// icosphere faceting gives neighbouring front points genuinely different
// per-point remaining distance d(p) (by design NOT smoothed -- only
// direction is smoothed, see the march step-length note above), so even
// a mild, uniform-radius contraction produces real per-edge length
// swings well above a "tight" bound like 1.5x; 1.5 measurably crushed
// the benign sphere case's average step to ~35% of its own d(p) and
// left residualMax close to a full t (should be near-zero for a benign
// field per the design doc). Loosened so the clamp only catches
// genuinely severe local distortion (sharp features), not ordinary
// faceting noise.
constexpr double kStretchMax = 6.0;
// Land exactly on the STL once the post-clamp remaining distance is
// under this fraction of t. MEASURED (offset_sphere): 0.05 left
// residualMax slightly above the gate's 0.1*t bar (icosphere faceting
// noise routinely leaves ~5-15% of t unmarched even on a benign
// field); 0.2 lands those points exactly instead of leaving them just
// short, which is exactly what "landing" is for.
constexpr double kLandFrac = 0.2;
constexpr double kMinHeightFrac = 0.05;    // min prism height, fraction of t, else invalid
// ACHIEVED-vs-REQUESTED layer height floor: a prism whose achieved mean
// height (volume / top area) falls below this fraction of the layer
// step the dict ASKED FOR is refused, its stack condemned, and the wall
// keeps the offset cut. See prismValid for the defect that forced it
// (bm_solver_wfp_interfoam's 0.88 mm first layer against a requested
// 3.1 mm) and for why the three pre-existing floors cannot see it.
//
// THE THRESHOLD IS MEASURED, not chosen. The per-prism
// achieved/requested ratio was dumped for every EMITTED wall-layer
// prism (NINJA_LAYER_HEIGHT_HIST) on all four real layered cases. Both
// shipped grids give the same shape: a dense population from 0.35
// upward (bm_layers_wfp alone has 38 323 prisms in the 0.35-0.40 bin,
// bm_layers_fishpassage_coarse 24 346 in 0.75-0.80) over a thin,
// flat tail below it (~500-700 prisms per 0.05-wide bin on wfp, ~50-400
// on fishpassage). Cumulative cost of a floor, as a fraction of all
// first-layer prisms:
//     floor   fishpassage_coarse   wfp
//     0.30          0.21 %         0.38 %
//     0.35          0.30 %         0.84 %
//     0.45          0.35 %         7.1 %      <- eats the dense band
// The located defect (bm_solver_wfp_interfoam cell 1052979) sits at
// 0.876 / 3.1 = 0.283. 0.35 is the LARGEST floor that still excludes
// the dense band entirely, and it clears the defect by 24 % rather
// than by the 6 % a 0.30 floor would leave. Under 1 % of first-layer
// prisms on every shipped grid, which is the bar this had to meet.
constexpr double kMinAchievedHeightFrac = 0.35;
// Measurement/tuning override for the constant above
// (NINJA_MIN_ACHIEVED_HEIGHT_FRAC). Opt-in, byte-identical when unset;
// it exists so the threshold sweep below can be re-run on a future
// case without a rebuild.
double minAchievedHeightFrac() {
    static const double v = [] {
        const char* e = std::getenv("NINJA_MIN_ACHIEVED_HEIGHT_FRAC");
        if (!e || !*e) return kMinAchievedHeightFrac;
        return std::atof(e);
    }();
    return v;
}
// Opt-in per-prism achieved/requested height histogram
// (NINJA_LAYER_HEIGHT_HIST), zero cost when unset. Bins 0..19 are
// achieved/requested in [0,1) at 0.05 apiece; bin 20 is >= 1. Rank-0
// serial accumulator (the validation loop is not threaded), printed by
// applyLayers.
// [layer index j (innermost = 0), capped at 9][ratio bin]
std::array<std::array<long, 21>, 10> g_heightHistBins{};
bool layerHeightHist() {
    static const bool on = [] {
        const char* e = std::getenv("NINJA_LAYER_HEIGHT_HIST");
        return e && *e && std::string(e) != "0";
    }();
    return on;
}
// checkMesh's own high-aspect-ratio threshold
// (primitiveMesh::aspectThreshold_ = 1000). A prism above it is
// refused by prismValid, so the march never emits a cell checkMesh
// would call an error. Dimensionless, not an absolute epsilon.
constexpr double kMaxPrismAspect = 1000.0;
// prismValid's bottom-fan test: a fan triangle counts as reflex only when
// its NEGATIVE area exceeds this fraction of the face's. Dimensionless.
constexpr double kReflexAreaFrac = 0.01;
// Micro-edge rigid clusters (measured on offset_sphere + merged
// slivers): front points joined by an edge shorter than this
// fraction of the current layer step march RIGIDLY (shared direction,
// shared step length). Rationale: across a micro edge (e.g. a merged
// sliver's ~1e-3 cut facet vs a ~2e-2 step) ANY differential march is
// a huge relative distortion -- the stretch clamp crushes those points
// to ~0 while their neighbours march in full, the lopsided quads fail
// the twist test (measured: 63/152 seed failures, all at or adjacent
// to micro facets). Marching the cluster as one rigid body removes the
// distortion by construction: micro facets extrude as thin-but-valid
// prisms and their neighbours see no relative motion across the shared
// edge.
constexpr double kMicroEdgeFrac = 0.25;
constexpr int kNonInversionMaxIter = 10;   // bisection outer-iteration cap
constexpr double kFreezeScaleThreshold = 0.05; // scaleFactor below this counts as "frozen"

Vec3 normalizeOrZero(const Vec3& v) {
    const double n = norm(v);
    if (n < 1e-300) {
        return Vec3{0, 0, 0};
    }
    return v * (1.0 / n);
}

Vec3 newellNormal(const std::vector<Vec3>& loop) {
    Vec3 nw{0, 0, 0};
    const std::size_t n = loop.size();
    for (std::size_t i = 0; i < n; ++i) {
        const Vec3& pc = loop[i];
        const Vec3& pn = loop[(i + 1) % n];
        nw.x += (pc.y - pn.y) * (pc.z + pn.z);
        nw.y += (pc.z - pn.z) * (pc.x + pn.x);
        nw.z += (pc.x - pn.x) * (pc.y + pn.y);
    }
    return nw;
}

double loopArea(const std::vector<Vec3>& loop) { return 0.5 * norm(newellNormal(loop)); }

Vec3 loopCentroid(const std::vector<Vec3>& loop) {
    Vec3 c{0, 0, 0};
    for (const Vec3& p : loop) c = c + p;
    return c * (1.0 / static_cast<double>(loop.size()));
}

// Sub-volume (tet-fan-about-centroid) validity check for one prism:
// top loop (outward-from-core orientation, as stored on the mesh) +
// bottom loop (same point order, marched positions). Positive
// sub-volumes for EVERY fan triangle of every face, oriented outward
// from the prism's own centroid, plus a minimum-height floor.
// failCode (debug instrumentation, optional): 0 = valid, 1 = degenerate
// loop, 2 = minEdge below minHeight, 3 = side-quad twist, 4 = bottom
// non-convex fan, 5 = non-positive prism volume, 6 = bottom face folded
// (normal not pointing away from the top), 7 = side quad pyramid
// inverted (outward normal points into the prism), 8 = mean height
// (volume / top area) below minHeight (sheared-flat prism), 9 = cell
// aspect ratio above maxAspect (checkMesh's own metric; see below),
// 10 = ACHIEVED mean height below minAchievedFrac of the REQUESTED
// step thickness (collapsed prism; see below).
//
// `meanHeightOut` (optional) receives the achieved mean height
// (volume / top area) whenever the prism reaches the height checks --
// the measurement the achieved-vs-requested guard and its gate key are
// built on. Left untouched when an earlier check rejects the prism.
bool prismValid(const std::vector<Vec3>& top, const std::vector<Vec3>& bot, double minHeight, double minVolEps,
                int* failCode = nullptr, double maxAspect = 0.0, double minAchievedHeight = 0.0,
                double* meanHeightOut = nullptr) {
    if (failCode) *failCode = 0;
    const std::size_t n = top.size();
    if (n != bot.size() || n < 3) {
        if (failCode) *failCode = 1;
        return false;
    }
    double minEdge = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < n; ++i) {
        minEdge = std::min(minEdge, norm(bot[i] - top[i]));
    }
    if (minEdge < minHeight) {
        if (failCode) *failCode = 2;
        return false;
    }

    // Twist/planarity guard (MEASURED necessary, offset_sphere):
    // checkMesh's own face-pyramid orientation check triangulates each
    // face independently (a stricter, per-face test than this
    // function's own aggregate prism-volume check below) and flagged
    // ~140 "wrongly oriented" internal/boundary faces on the first
    // working march -- all on the march's OWN side quads, none on the
    // untouched top faces. Root cause: a quad whose two fan triangles
    // (from a shared anchor corner) disagree in normal direction is
    // legitimately ambiguous ("twisted") for a per-triangle consumer
    // like checkMesh's pyramid decomposition, even though this
    // function's own aggregate prism-volume check (below) doesn't care.
    // MEASURED: two looser alternatives (comparing the two DIAGONALS'
    // sums; anchoring against the quad's own Newell normal) both let
    // through geometry that still failed checkMesh -- this strict
    // same-anchor-pair test is the one that measurably produced a
    // clean "Mesh OK" + integrity pass. Disclosed cost: it is
    // conservative enough to also reject some quads a human would call
    // "fine", so it (via the non-inversion bisection loop below)
    // shrinks the march noticeably short of a full landing on this
    // benchmark's icosphere faceting -- reported, not silently
    // tolerated: residualMax is measured above the gate's target on
    // this case and reported as such.
    for (std::size_t i = 0; i < n; ++i) {
        const Vec3& p0 = top[i];
        const Vec3& p1 = top[(i + 1) % n];
        const Vec3& p2 = bot[(i + 1) % n];
        const Vec3& p3 = bot[i];
        const Vec3 nA = cross(p1 - p0, p2 - p0);
        const Vec3 nB = cross(p2 - p0, p3 - p0);
        if (dot(nA, nB) < 0.0) {
            if (failCode) *failCode = 3;
            return false;
        }
    }
    {
        // Reference = the loop's own Newell normal (the first fan triangle is
        // pure noise when bot[1] is a collinear hanging node -- every wall
        // face at a refinement-level transition has one), and a fan triangle
        // only counts as reflex when its negative area is a meaningful
        // fraction of the face (relative, not an exact sign test). The
        // non-inversion clamp uses the same test: with the exact one it
        // halved every point of such a face on a coin flip and crushed the
        // healthy prisms sharing them (MEASURED, hydrofoil window: 828 of 901
        // collapsed stacks, 1831 -> 799 dropped faces).
        const Vec3 bN = newellNormal(bot); // |bN| ~ 2 * area
        const double bN2 = dot(bN, bN);
        for (std::size_t i = 1; i + 1 < n; ++i) {
            const Vec3 bni = cross(bot[i] - bot[0], bot[(i + 1) % n] - bot[0]);
            if (dot(bN, bni) < -kReflexAreaFrac * bN2) {
                if (failCode) *failCode = 4;
                return false;
            }
        }
    }
    // Absolute bottom-face orientation (measured necessary on
    // offset_rotcube): at a convex STL edge the landing snap can fold a whole
    // bottom loop past the edge plane, flipping its normal INTO the prism.
    // The fan check above only tests internal consistency and the
    // aggregate volume below survives a small local flip, so checkMesh's
    // per-face pyramid test was the first thing to notice (3 misoriented
    // wall faces on the rotcube's rounded bottom edge). The bottom's
    // outward normal must point away from the top (the march direction).
    {
        const Vec3 botNw = newellNormal(bot);
        const Vec3 topC = loopCentroid(top);
        const Vec3 botC = loopCentroid(bot);
        if (dot(botNw, botC - topC) <= 0.0) {
            if (failCode) *failCode = 6;
            return false;
        }
    }

    Vec3 centroid{0, 0, 0};
    for (const Vec3& p : top) centroid = centroid + p;
    for (const Vec3& p : bot) centroid = centroid + p;
    centroid = centroid * (1.0 / static_cast<double>(2 * n));

    auto tetVol = [&](const Vec3& a, const Vec3& b, const Vec3& c) {
        return dot(a - centroid, cross(b - centroid, c - centroid)) / 6.0;
    };

    // Sub-volumes summed (divergence theorem, anchored at the prism's
    // own centroid): a genuine tet-fan-about-centroid decomposition,
    // exactly as the brief specifies. NOTE (disclosed simplification):
    // a per-face fan-triangulation choice (which diagonal of a quad, or
    // which vertex a polygon fans from) is not unique for a NON-PLANAR
    // face -- and both the offset cut's own top faces and this march's
    // side quads (top and bottom corners moving by different amounts)
    // are routinely non-planar. Requiring every INDIVIDUAL fan
    // triangle's sub-tet to be positive (the first implementation here)
    // measurably rejected plainly valid, merely-twisted prisms as
    // "inverted" purely from that triangulation-anchor artifact
    // (MEASURED: >90% of a benign sphere front spuriously failing).
    // The robust, anchor-invariant criterion is the SUM of every sub-
    // tet (equivalently, the standard mesh-volume-via-divergence-
    // theorem formula applied to the prism's own closed, consistently-
    // outward-oriented face set) -- this is a true measure of the
    // prism's actual (possibly warped-but-still-simple) volume, and
    // still catches genuine self-intersection/inversion (which drives
    // the total negative), just not spurious triangulation artifacts.
    double totalVol = 0.0;
    auto fanFromCentroid = [&](const std::vector<Vec3>& loop, bool reverseOrientation) {
        const Vec3 fc = loopCentroid(loop);
        const std::size_t m = loop.size();
        for (std::size_t i = 0; i < m; ++i) {
            Vec3 a = fc;
            Vec3 b = loop[i];
            Vec3 cc = loop[(i + 1) % m];
            if (reverseOrientation) std::swap(b, cc);
            totalVol += tetVol(a, b, cc);
        }
    };
    // Top face: outward-FROM-PRISM is the reverse of outward-from-core
    // (the stored loop direction).
    fanFromCentroid(top, true);
    // Bottom face: outward-from-prism == the stored top-loop direction
    // (bottom is roughly parallel, further from the core).
    fanFromCentroid(bot, false);
    // Side quads in their TOPOLOGICAL outward orientation
    // {top[i], top[i+1], bot[i+1], bot[i]} -- the order the emitter
    // writes (see the side-face emission below). MEASURED (bm_layers_wfp,
    // 31 open cells / 44 misoriented pyramids that check_integrity's
    // undirected-edge closure could not see): this used to resolve each
    // quad's winding GEOMETRICALLY (flip if its Newell normal pointed at
    // the centroid), which silently "repaired" a folded quad here, in
    // the quality evaluator and in the emitter alike -- so an inverted
    // side face passed every guard and was written with a winding
    // inconsistent with its neighbours, i.e. an open cell. Orientation
    // is a property of the topology, never of the geometry; a quad
    // whose outward normal points inward is a defect to be caught, not
    // re-wound.
    for (std::size_t i = 0; i < n; ++i) {
        const std::vector<Vec3> quad{top[i], top[(i + 1) % n], bot[(i + 1) % n], bot[i]};
        totalVol += tetVol(quad[0], quad[1], quad[2]);
        totalVol += tetVol(quad[0], quad[2], quad[3]);
        // Per-face pyramid sign (checkMesh's face-pyramid test, anchor-
        // invariant: Newell area vector dotted with centroid-to-face-
        // centre) -- a folded side quad must fail here even when the
        // aggregate volume below still comes out positive.
        if (dot(newellNormal(quad), loopCentroid(quad) - centroid) <= 0.0) {
            if (failCode) *failCode = 7;
            return false;
        }
    }
    if (totalVol <= minVolEps) {
        if (failCode) *failCode = 5;
        return false;
    }
    // Mean-height floor: volume over top-face area. The per-point
    // |bot - top| floor above is blind to SHEAR -- MEASURED
    // (bm_layers_wfp): prisms whose bottom points all moved >= minHeight
    // but almost tangentially, leaving a cell 1e-6..3e-5 m thick under a
    // 0.05-0.2 m face (checkMesh aspect ratio 5933, the 1e-7 m^2 seam
    // quads behind its max skewness 221). Anchor-invariant, no new
    // tolerance: the same minHeight, applied to the mean height.
    {
        const double topArea = 0.5 * norm(newellNormal(top));
        const double meanHeight = topArea > 1e-300 ? totalVol / topArea : 0.0;
        if (meanHeightOut && topArea > 1e-300) *meanHeightOut = meanHeight;
        if (minHeight > 0.0 && topArea > 1e-300 && meanHeight < minHeight) {
            if (failCode) *failCode = 8;
            return false;
        }
        // ACHIEVED-vs-REQUESTED floor (the collapsed-prism guard).
        //
        // MEASURED (bm_solver_wfp_interfoam, f6e3d05): cell 1052979 at
        // (3.546 7.017 5.355) is a prism 0.0363 x 0.0386 x 0.000876 m,
        // V = 1.68e-07 m3 -- 0.88 mm of layer where the dict's 6-layer
        // r=1.8 stack asks for a 3.1 mm first layer, i.e. 28% of what
        // was requested. checkMesh calls it "Cell volumes OK" and "Max
        // skewness 3.98 OK"; interFoam does not: Co = 0.5*dt*sum|phi|/V
        // over 31 sub-1e-06 m3 cells forced the controller to cut
        // deltaT 20x (2.8e-03 -> 1e-04) and still reported max Courant
        // 2.79 against the case's own maxCo 2. A prism squeezed to a
        // quarter of its requested height is a bad cell that happens to
        // pass checkMesh, and the contract is that the mesh never ships
        // one: where the layer does not fit, drop the stack and stay at
        // the cut.
        //
        // WHY THE OTHER THREE FLOORS MISS IT, each for its own reason:
        //   * the per-point |bot - top| floor and the mean-height floor
        //     above are both kMinHeightFrac = 5% of the requested step,
        //     and 28% clears 5% comfortably;
        //   * the aspect ceiling (aadfd00) measures height against
        //     WIDTH: 0.876 mm under a 36 x 39 mm face is aspect ~44,
        //     nowhere near checkMesh's 1000.
        // This test is the remaining axis and the one the defect is
        // actually on -- achieved height against the height that was
        // ASKED FOR. It is a dimensionless ratio of two lengths that
        // scale together (both are the dict's own thickness schedule),
        // so unlike aadfd00's absolute-length floors it carries its
        // context with it: coarsening the grid does not move it.
        if (minAchievedHeight > 0.0 && topArea > 1e-300 && meanHeight < minAchievedHeight) {
            if (failCode) *failCode = 10;
            return false;
        }
    }
    // Aspect-ratio ceiling, checkMesh's OWN metric, transcribed from
    // primitiveMeshTools::cellClosedness rather than paraphrased:
    //     sumMagClosed[dir] = sum over the cell's faces of |Sf[dir]|
    //     aspect = max( max(sumMagClosed)/min(sumMagClosed),
    //                   (1/6)*cmptSum(sumMagClosed) / V^(2/3) )
    // and `maxAspect` is checkMesh's own default threshold.
    //
    // WHY, and why the two floors above do not cover it (MEASURED,
    // bm_layers_descargador at grid factor 8): both minHeight floors are
    // relative to the requested layer THICKNESS (kMinHeightFrac * tStep)
    // and blind to how WIDE the wall face is. Nothing ties the two
    // together -- a stack can be thin against the face it grows on
    // whenever the thickness is held while the grid coarsens (which
    // coarsen_sweep.sh does deliberately, and which the dict could ask
    // for directly), or simply because expansionRatio drives the wall
    // layer far below the rest. Measured: a 0.09 m stack on a 1 m wall
    // face gives a prism 8e-4 m tall that passes the height floor
    // (8e-4 > 0.05 * 0.011) and lands at aspect ratio 1258.
    // checkMesh reported 5 such cells at max 1942.68 -- an ERROR, i.e.
    // a mesh with a quality problem, which is exactly what asking for
    // cells too big for the requested layers must never produce.
    //
    // Refusing is the robust direction and the one the rest of this
    // stage already takes: the prism is reverted at this step, its stack
    // is condemned, and the wall keeps the offset cut -- fewer layers,
    // mesh away from the wall, never an out-of-floor cell.
    if (maxAspect > 0.0) {
        Vec3 sumMag{0, 0, 0};
        auto addFaceMag = [&](const std::vector<Vec3>& loop) {
            const Vec3 sf = newellNormal(loop) * 0.5;
            sumMag = sumMag + Vec3{std::fabs(sf.x), std::fabs(sf.y), std::fabs(sf.z)};
        };
        addFaceMag(top);
        addFaceMag(bot);
        for (std::size_t i = 0; i < n; ++i) {
            addFaceMag(std::vector<Vec3>{top[i], top[(i + 1) % n], bot[(i + 1) % n], bot[i]});
        }
        const double mx = std::max(sumMag.x, std::max(sumMag.y, sumMag.z));
        const double mn = std::min(sumMag.x, std::min(sumMag.y, sumMag.z));
        // OpenFOAM divides by (min + ROOTVSMALL); the same guard, so a
        // zero-projection direction gives a huge ratio rather than a NaN.
        double aspect = mx / (mn + 1e-150);
        const double cmptSum = sumMag.x + sumMag.y + sumMag.z;
        aspect = std::max(aspect, (1.0 / 6.0) * cmptSum / std::cbrt(totalVol * totalVol));
        if (aspect > maxAspect) {
            if (failCode) *failCode = 9;
            return false;
        }
    }
    return true;
}

// --- Opt-in layer-quality instrumentation ----------------------------
// (NINJA_LAYER_QUALITY_STATS, opt-in, zero cost when
// unset). Reimplements, locally, the two
// checkMesh metrics the march's known defect family fails
// (primitiveMeshCheck semantics): face-pyramid orientation and
// centroid-decomposition tet quality. Both share the SAME tet set --
// per face, fan the face's own loop from its centroid, each fan
// triangle forms a tet with the CELL's centroid as apex (exactly
// prismValid's own centroid-tet decomposition, generalized to report
// PER-FACE and PER-TET results instead of prismValid's aggregate
// sum). Pyramid orientation = sign of a face's tet-sum (checkMesh's
// owner/neighbour pyramid test, cell-local: correct orientation from
// a cell's own point of view does not depend on whether the face is
// this cell's boundary, or an internal face's owner or neighbour
// side); tet quality = sign of each INDIVIDUAL tet (checkMesh's
// stricter per-tet centroid-decomposition test -- this is the one
// prismValid's own comment above explains was rejected as the
// PRISM-validity criterion because it flags plainly-valid twisted
// quads, but that is exactly the failure mode this diagnostic exists
// to surface, not hide).
struct LayerCellQuality {
    int nBadPyramidFaces = 0;
    int nBadTetFaces = 0;
    int nHighSkewFaces = 0; // checkMesh faceSkewness above its own thresholds (see evaluateLayerCellQuality)
    int nHighNonOrthFaces = 0; // non-orthogonality above kNonOrthMaxDeg (see evaluateLayerCellQuality)
    double worstTetNorm = std::numeric_limits<double>::max(); // min tetPointRef::quality() over all tets (dimensionless)
    Vec3 centroid{0, 0, 0}; // this prism's own centroid (cached for the NEXT step's owner-side check)
};

// OpenFOAM's tetPointRef::mag()/circumRadius()/quality()
// (src/OpenFOAM/meshes/primitiveShapes/tetrahedron/tetrahedronI.H),
// transcribed exactly, vertex order (a,b,c,d). quality() is the
// (signed) volume normalized by the volume of the circumscribing
// sphere's inscribed regular tet -- scale-invariant, quality 1 for a
// perfect regular tet, near 0 or negative for degenerate/inverted
// ones. polyMeshTetDecomposition::minTetQuality = sqr(SMALL).
double tetMag(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    return dot(b - a, cross(c - a, d - a)) / 6.0;
}
double tetCircumRadius(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    const double GREAT = 1e15;
    const double ROOTVSMALL = 1e-150;
    const Vec3 va = b - a;
    const Vec3 vb = c - a;
    const Vec3 vc = d - a;
    const double lambda = dot(vc, vc) - dot(va, vc);
    const double mu = dot(vb, vb) - dot(va, vb);
    const Vec3 ba = cross(vb, va);
    const Vec3 ca = cross(vc, va);
    const Vec3 num = ba * lambda - ca * mu;
    const double denom = dot(vc, ba);
    if (std::fabs(denom) < ROOTVSMALL) return GREAT;
    return norm(va * 0.5 + num * (0.5 / denom));
}
double tetQuality(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    const double GREAT = 1e15;
    const double ROOTVSMALL = 1e-150;
    const double R = std::min(tetCircumRadius(a, b, c, d), GREAT);
    return tetMag(a, b, c, d) / (8.0 / (9.0 * std::sqrt(3.0)) * R * R * R + ROOTVSMALL);
}
// OpenFOAM: polyMeshTetDecomposition::minTetQuality = sqr(SMALL).
const double kMinTetQuality = 1e-30;

// polyMeshTetDecomposition::minQuality (polyMeshTetDecomposition.C):
// fans a face's OTHER points (walking the loop starting at
// `basePtIdx`, skipping it) into tets sharing the edge (cellCentre,
// loop[basePtIdx]); returns the WORST (min) tetQuality() over that
// fan. `isOwner` mirrors OpenFOAM's own/nei point-order swap (own:
// (facePt, nextFacePt); nei: (nextFacePt, facePt)).
double baseFanMinQuality(const Vec3& cC, const std::vector<Vec3>& loop, std::size_t basePtIdx, bool isOwner) {
    const std::size_t m = loop.size();
    const Vec3& tetBasePt = loop[basePtIdx];
    double worst = std::numeric_limits<double>::max();
    for (std::size_t tetPtI = 1; tetPtI + 1 < m; ++tetPtI) {
        const std::size_t facePtI = (tetPtI + basePtIdx) % m;
        const std::size_t otherFacePtI = (facePtI + 1) % m;
        const Vec3& pFace = loop[facePtI];
        const Vec3& pOther = loop[otherFacePtI];
        const Vec3& pA = isOwner ? pFace : pOther;
        const Vec3& pB = isOwner ? pOther : pFace;
        worst = std::min(worst, tetQuality(cC, tetBasePt, pA, pB));
    }
    return worst;
}

// polyMeshTetDecomposition::findSharedBasePoint / findBasePoint: TRUE
// iff SOME face vertex, used as the fan's shared base point, gives a
// fan (both sides, if `nCcPtr` given -- internal face; owner side
// only, if `nCcPtr` is null -- boundary face) whose worst tet quality
// clears `tol` -- i.e. a usable tet decomposition exists for this face
// AT ALL. This (not the simple centroid fan above) is checkMesh's own
// dominant discriminator in practice: a centroid-anchored fan stays
// the right sign even for a merely-warped-but-simple face (its
// signed-volume magnitude just gets small), while a fan anchored on
// one of the face's OWN corners can go genuinely inverted once the
// polygon is non-planar/concave enough -- exactly the defect family
// this diagnostic targets, and MEASURED (see report) to be the check
// that actually reproduces checkMesh's lowQualityTetFaces counts,
// where the plain centroid-fan sign test does not.
bool faceHasUsableBasePoint(const Vec3& oCc, const Vec3* nCcPtr, const std::vector<Vec3>& loop, double tol) {
    const std::size_t m = loop.size();
    for (std::size_t bp = 0; bp < m; ++bp) {
        const double ownQ = baseFanMinQuality(oCc, loop, bp, true);
        if (!nCcPtr) {
            if (ownQ > tol) return true;
            continue;
        }
        const double neiQ = baseFanMinQuality(*nCcPtr, loop, bp, false);
        if (std::min(ownQ, neiQ) > tol) return true;
    }
    return false;
}

// Area-weighted polygon centroid (fan-triangulated from the plain
// vertex average): needed for a faithful cell-centroid pyramid
// decomposition (loopCentroid's unweighted vertex average is fine for
// near-regular quads but biases badly on the warped, irregular
// polygons this diagnostic exists to catch).
Vec3 areaWeightedCentroid(const std::vector<Vec3>& loop) {
    const Vec3 pAvg = loopCentroid(loop);
    const std::size_t m = loop.size();
    Vec3 sumAreaC{0, 0, 0};
    double sumArea = 0.0;
    for (std::size_t i = 0; i < m; ++i) {
        const Vec3& a = loop[i];
        const Vec3& b = loop[(i + 1) % m];
        const double area = 0.5 * norm(cross(a - pAvg, b - pAvg));
        const Vec3 triC = (pAvg + a + b) * (1.0 / 3.0);
        sumAreaC = sumAreaC + triC * area;
        sumArea += area;
    }
    return sumArea > 1e-300 ? sumAreaC * (1.0 / sumArea) : pAvg;
}

// One face's centroid-decomposition tet fan, evaluated against an
// explicit APEX (a cell centroid, owner or neighbour) -- shared by the
// per-prism evaluator below and the extra owner-side (core cell) top-
// face check, since checkMesh's pyramid/tet tests are the SAME
// arithmetic from either side of a face, just with a different apex
// and (for the non-owner side) a reversed loop.
// Per checkFaceTets (polyMeshTetDecomposition.C): for each consecutive
// face-point pair (a,b) the tet is (a, b, faceCentre, cellCentre) --
// i.e. tetPointRef(P0=apex, P1=fc, P2=loop[i], P3=loop[i+1]) once the
// scalar-triple-product identity a.(bxc) = c.(axb) is used to move the
// apex into the P0 slot; this loop's own orientation must already be
// "outward from apex's cell" (checkMesh always uses the mesh's stored,
// owner-outward face point order and varies only which cell's centre
// is used as apex -- an owner-side and a reversed neighbour-side tet
// fan are algebraically equivalent to that, see report). tetQuality()
// is checkMesh's actual metric (not a raw normalized volume); a tet is
// bad when quality <= minTetQuality, matching checkFaceTets' own
// (owner: mag > -tol) / (neighbour: mag < tol) pair once the sign
// flip from the apex-in-P0 reindexing is folded in.
// Returns true iff this apex's centroid fan alone already condemns the
// face (kept as a bool return, not an r.nBadTetFaces++, so the caller
// can OR it together with the Part-3 base-point search below and
// increment the per-face tet counter exactly once per face).
bool evalFaceAgainstApex(const std::vector<Vec3>& loop, const Vec3& apex, LayerCellQuality& r) {
    const Vec3 fc = areaWeightedCentroid(loop);
    const std::size_t m = loop.size();
    double faceSum = 0.0;
    bool faceBadTet = false;
    for (std::size_t i = 0; i < m; ++i) {
        const double q = tetQuality(apex, fc, loop[i], loop[(i + 1) % m]);
        faceSum += tetMag(apex, fc, loop[i], loop[(i + 1) % m]);
        r.worstTetNorm = std::min(r.worstTetNorm, q);
        if (q <= kMinTetQuality) faceBadTet = true;
    }
    if (faceSum <= 0.0) ++r.nBadPyramidFaces;
    return faceBadTet;
}

// checkMesh's face skewness (primitiveMeshTools::faceSkewness /
// boundaryFaceSkewness, transcribed): the in-plane offset of the face
// centre from the owner->neighbour line (internal) or from the owner
// centre's normal projection (boundary), normalised by the face's own
// half-extent in that direction. `loop` outward from the owner
// (`ownCc`); `neiCc` null for the boundary form. MEASURED necessary
// (bm_layers_wfp): a prism whose bottom points moved mostly
// tangentially keeps a positive volume and valid pyramids, but its
// SEAM side quads lie almost flat against the wall -- the prism
// centroid projects ~w/2 outside a face only h thick, boundary
// skewness 200+ (checkMesh's max skewness 221 on that case).
double faceSkewnessOf(const std::vector<Vec3>& loop, const Vec3& ownCc, const Vec3* neiCc) {
    const double ROOTVSMALL = 1e-150;
    const Vec3 fc = areaWeightedCentroid(loop);
    const Vec3 fA = newellNormal(loop) * 0.5;
    const Vec3 Cpf = fc - ownCc;
    Vec3 sv;
    double fd;
    if (neiCc) {
        const Vec3 d = *neiCc - ownCc;
        sv = Cpf - d * (dot(fA, Cpf) / (dot(fA, d) + ROOTVSMALL));
        fd = 0.2 * norm(d) + ROOTVSMALL;
    } else {
        const Vec3 n = normalizeOrZero(fA);
        const Vec3 d = n * dot(n, Cpf);
        sv = Cpf - d;
        fd = 0.4 * norm(d) + ROOTVSMALL;
    }
    const Vec3 svHat = sv * (1.0 / (norm(sv) + ROOTVSMALL));
    for (const Vec3& p : loop) fd = std::max(fd, std::fabs(dot(svHat, p - fc)));
    return norm(sv) / fd;
}
// checkMesh's own threshold (primitiveMesh::skewThreshold_ = 4, applied
// to internal and boundary faces alike) -- the consumer's number, not a
// tuned knob of this file.
constexpr double kSkewMax = 4.0;

// Non-orthogonality of a prism face, in degrees: the angle between the face
// normal and the line from `ownCc` to the neighbour centre (`neiCc`) or, where
// the neighbour is not decided yet, to the face centre -- the same one-sided
// proxy the skewness test above uses for side quads.
//
// No validity test sees it: a prism can be closed, positively oriented,
// unskewed and still share a side with a sibling whose centre lies ALONG that
// side. MEASURED, win_wfp_trans (shipped resolution): where a diagonal wall
// meets the grid-aligned ceiling the wall faces are a zigzag chain of needle
// triangles, 140 x 6 mm; their prisms meet on the short sides with centres
// displaced along the strip, 17 faces above 70 degrees, worst 82.8 (the proxy
// gives 81 on that face; a square prism gives 0, a right-triangle one 26).
// Layers off, the same window tops out at 52.8.
double faceNonOrthDegOf(const std::vector<Vec3>& loop, const Vec3& ownCc, const Vec3* neiCc) {
    const Vec3 n = newellNormal(loop);
    const Vec3 d = (neiCc ? *neiCc : areaWeightedCentroid(loop)) - ownCc;
    const double den = norm(n) * norm(d);
    if (den < 1e-300) return 0.0;
    const double c = std::max(-1.0, std::min(1.0, dot(n, d) / den));
    return std::acos(c) * (180.0 / 3.14159265358979323846);
}
// checkMesh's "severely non-orthogonal" threshold (primitiveMesh::nonOrthThreshold_ = 70).
constexpr double kNonOrthMaxDeg = 70.0;
// Deepest a landed point may sit inside the solid and still be emitted
// (see the seal guard), in units of the face's h_local. Measured maximum
// 0.2 h on a hydrofoil's sharp edges; the smoothed field bounds it.
constexpr double kMaxBuriedDepthFrac = 0.25;

// Volume-weighted polyhedron centroid (same pass1/pass2 scheme as
// cellCentroidOf below), given the cell's faces as loops ALREADY
// oriented outward from the cell -- the case for a freshly-built
// prism, where the outward orientation of each face is known directly
// from the march's own bookkeeping instead of an owner/neighbour flag.
// MEASURED necessary (see report): OpenFOAM cell centres
// (primitiveMesh::cellCentres) are this volume-weighted centroid, not
// a vertex average, and the two diverge enough on the warped prisms
// this diagnostic targets to change tet-quality verdicts.
Vec3 polyhedronCentroidFromOutwardLoops(const std::vector<std::vector<Vec3>>& loops) {
    Vec3 est{0, 0, 0};
    int nf = 0;
    for (const auto& loop : loops) {
        est = est + areaWeightedCentroid(loop);
        ++nf;
    }
    if (nf == 0) return est;
    est = est * (1.0 / static_cast<double>(nf));

    // primitiveMesh::makeCellCentresAndVols, transcribed: per face, a
    // pyramid of 3*volume = fA . (fc - est) and centroid 3/4 fc + 1/4 est
    // (NOT a per-vertex tet fan -- the two differ on non-planar faces,
    // which the thin sheared prisms this evaluator must see are).
    double sumVol = 0.0;
    Vec3 sumVolC{0, 0, 0};
    for (const auto& loop : loops) {
        const Vec3 fc = areaWeightedCentroid(loop);
        const Vec3 fA = newellNormal(loop) * 0.5;
        const double pyr3Vol = dot(fA, fc - est);
        const Vec3 pc = fc * 0.75 + est * 0.25;
        sumVol += pyr3Vol;
        sumVolC = sumVolC + pc * pyr3Vol;
    }
    return std::fabs(sumVol) > 1e-300 ? sumVolC * (1.0 / sumVol) : est;
}

// Evaluates a newly-created prism cell's OWN faces (top from the
// prism/neighbour side, bottom, sides) against its own centroid.
// `ownerCentroid` (null when there is no such cell): the pre-existing
// "core" cell across the top face, i.e. the OTHER side of that internal
// face. MEASURED necessary: the known defect family lives on the top
// face's OWNER-side pyramid (the core cell, which can be a large or
// merged/non-convex cut cell whose own centroid sits far enough from
// a warped wall pentagon to flip that one face's sign), which the
// prism's own well-shaped local centroid alone never reproduces.
LayerCellQuality evaluateLayerCellQuality(const std::vector<Vec3>& top, const std::vector<Vec3>& bot,
                                           double /*lengthScale*/, const Vec3* ownerCentroid,
                                           const std::vector<const Vec3*>* sideNeighbourCentroids = nullptr) {
    LayerCellQuality r;
    const std::size_t n = top.size();
    if (n != bot.size() || n < 3) return r;

    // Top face, NEIGHBOUR (prism) side: outward-FROM-PRISM is the
    // reverse of the stored (outward-from-core) loop direction.
    std::vector<Vec3> topRev(top.rbegin(), top.rend());
    // Bottom face: outward-from-prism == the stored top-loop direction.
    // Side quads in their topological outward orientation, exactly as
    // emitted (see prismValid's note: winding is never resolved
    // geometrically, so a folded quad is SEEN here as a bad pyramid).
    std::vector<std::vector<Vec3>> sideQuads(n);
    for (std::size_t i = 0; i < n; ++i) {
        sideQuads[i] = std::vector<Vec3>{top[i], top[(i + 1) % n], bot[(i + 1) % n], bot[i]};
    }

    // OpenFOAM cell centres (primitiveMesh::cellCentres) are the
    // volume-weighted centroid of the cell's own (outward-oriented)
    // faces -- MEASURED necessary (see report) that this diverges
    // enough from a plain vertex average to flip tet-quality verdicts
    // on the warped prisms this diagnostic targets.
    std::vector<std::vector<Vec3>> outwardLoops;
    outwardLoops.reserve(n + 2);
    outwardLoops.push_back(topRev);
    outwardLoops.push_back(bot);
    for (const auto& q : sideQuads) outwardLoops.push_back(q);
    const Vec3 centroid = polyhedronCentroidFromOutwardLoops(outwardLoops);
    r.centroid = centroid;

    // Each face's verdict combines BOTH real checkFaceTets tests: the
    // centroid fan (Parts 1/2) OR'd with the base-point search (Part
    // 3, findSharedBasePoint/findBasePoint) -- either one alone
    // condemns the face, exactly as checkFaceTets does (it flags on
    // the first failing test it finds, never "un-flags"). The top face
    // is the only one here with a known SECOND (owner) apex, so it
    // alone gets the real two-sided base-point search (`top` = the
    // face's own, owner-oriented point order, unreversed -- OpenFOAM's
    // convention); every other face (bottom, sides) only has this
    // prism's own centroid tracked here, so gets the single-sided
    // (boundary-style) variant -- same documented single-apex scope as
    // the pyramid check above.
    bool topBad = evalFaceAgainstApex(topRev, centroid, r);
    // Top face, OWNER (core cell) side, if known: stored loop order IS
    // outward-from-core already.
    if (ownerCentroid) {
        topBad = evalFaceAgainstApex(top, *ownerCentroid, r) || topBad;
        topBad = topBad || !faceHasUsableBasePoint(*ownerCentroid, &centroid, top, kMinTetQuality);
    } else {
        topBad = topBad || !faceHasUsableBasePoint(centroid, nullptr, topRev, kMinTetQuality);
    }
    if (topBad) ++r.nBadTetFaces;

    bool botBad = evalFaceAgainstApex(bot, centroid, r);
    botBad = botBad || !faceHasUsableBasePoint(centroid, nullptr, bot, kMinTetQuality);
    if (botBad) ++r.nBadTetFaces;

    for (const auto& quad : sideQuads) {
        bool quadBad = evalFaceAgainstApex(quad, centroid, r);
        quadBad = quadBad || !faceHasUsableBasePoint(centroid, nullptr, quad, kMinTetQuality);
        if (quadBad) ++r.nBadTetFaces;
    }

    // Skewness. Top face: the internal form when the owner (core) is
    // known. Bottom face: the BOUNDARY form from this prism's own centroid.
    // Side quads: see sideNeighbour below.
    if (ownerCentroid) {
        if (faceSkewnessOf(top, *ownerCentroid, &centroid) > kSkewMax) ++r.nHighSkewFaces;
    } else if (faceSkewnessOf(topRev, centroid, nullptr) > kSkewMax) {
        ++r.nHighSkewFaces;
    }
    if (faceSkewnessOf(bot, centroid, nullptr) > kSkewMax) ++r.nHighSkewFaces;
    // A side quad whose neighbouring front face also carries a prism this
    // step is an INTERNAL face between the two, and checkMesh scores it
    // against both centres; only a rim quad keeps the one-sided form (a
    // neighbour reverted later is re-scored one-sided by the caller).
    // MEASURED: scoring every side quad one-sided dropped healthy stacks on
    // sheared fronts -- bm_layers_wfp 335 -> 296, hydrofoil window 799 -> 550.
    auto sideNeighbour = [&](std::size_t i) -> const Vec3* {
        return sideNeighbourCentroids && i < sideNeighbourCentroids->size() ? (*sideNeighbourCentroids)[i] : nullptr;
    };
    for (std::size_t i = 0; i < n; ++i) {
        if (faceSkewnessOf(sideQuads[i], centroid, sideNeighbour(i)) > kSkewMax) ++r.nHighSkewFaces;
    }
    // Non-orthogonality, same split (checkMesh scores internal faces only).
    if (ownerCentroid && faceNonOrthDegOf(top, *ownerCentroid, &centroid) > kNonOrthMaxDeg) ++r.nHighNonOrthFaces;
    for (std::size_t i = 0; i < n; ++i) {
        if (faceNonOrthDegOf(sideQuads[i], centroid, sideNeighbour(i)) > kNonOrthMaxDeg) ++r.nHighNonOrthFaces;
    }
    return r;
}

// Faithful (OpenFOAM primitiveMeshCentres-style) polyhedron centroid:
// pass 1 estimates a centre from face centroids, pass 2 accumulates
// signed tet-pyramid volumes/centroids about that estimate (face
// orientation resolved via owner/neighbour, matching the mesh's own
// outward convention), and returns the volume-weighted average. Used
// ONLY for the ORIGINAL (pre-layers) core cell at step 0 of each
// stack, where the cell can be an arbitrary cut/merged polyhedron, not
// a prism -- the cheap vertex-average used elsewhere in this file is
// not trustworthy there.
Vec3 cellCentroidOf(const GeneratedMesh& m, int c) {
    const IntSpan cfaces = m.cellFacesOf(c);
    Vec3 est{0, 0, 0};
    int nf = 0;
    for (int fi : cfaces) {
        const IntSpan fp = m.faces.pointsOf(fi);
        std::vector<Vec3> loop;
        loop.reserve(static_cast<std::size_t>(fp.size()));
        for (int p : fp) loop.push_back(m.points[static_cast<std::size_t>(p)]);
        est = est + areaWeightedCentroid(loop);
        ++nf;
    }
    if (nf == 0) return est;
    est = est * (1.0 / static_cast<double>(nf));

    double sumVol = 0.0;
    Vec3 sumVolC{0, 0, 0};
    for (int fi : cfaces) {
        const IntSpan fp = m.faces.pointsOf(fi);
        std::vector<Vec3> loop;
        loop.reserve(static_cast<std::size_t>(fp.size()));
        for (int p : fp) loop.push_back(m.points[static_cast<std::size_t>(p)]);
        const bool isOwner = m.faces.owner[static_cast<std::size_t>(fi)] == c;
        const Vec3 fc = areaWeightedCentroid(loop);
        // Same OpenFOAM pyramid formula as polyhedronCentroidFromOutwardLoops;
        // the stored winding is outward from the OWNER, so flip the sign
        // on the neighbour side.
        const Vec3 fA = newellNormal(loop) * (isOwner ? 0.5 : -0.5);
        const double pyr3Vol = dot(fA, fc - est);
        const Vec3 pc = fc * 0.75 + est * 0.25;
        sumVol += pyr3Vol;
        sumVolC = sumVolC + pc * pyr3Vol;
    }
    return std::fabs(sumVol) > 1e-300 ? sumVolC * (1.0 / sumVol) : est;
}

// Max vertex distance from the area-weighted (Newell-normal) best-fit
// plane of a face loop -- the wall seed-face non-planarity diagnostic
// column, reused for every flagged stack's original wall
// face.
double loopNonPlanarity(const std::vector<Vec3>& loop) {
    if (loop.size() < 3) return 0.0;
    const Vec3 c = loopCentroid(loop);
    const Vec3 nrm = normalizeOrZero(newellNormal(loop));
    double worst = 0.0;
    for (const Vec3& p : loop) worst = std::max(worst, std::fabs(dot(p - c, nrm)));
    return worst;
}

// Which domain planes (0..5: xmin,xmax,ymin,ymax,zmin,zmax) `p` lies on,
// to an absolute tolerance scaled by the domain's own extent (grid
// points are exact, so this only ever fires for genuine boundary
// coincidence, not numerical noise from the march itself, which is why
// it is evaluated once per point up front, before any marching).
std::vector<int> domainPlanesOf(const Vec3& p, const MeshConfig& cfg) {
    const double ex = std::max(1e-300, cfg.max.x - cfg.min.x);
    const double ey = std::max(1e-300, cfg.max.y - cfg.min.y);
    const double ez = std::max(1e-300, cfg.max.z - cfg.min.z);
    const double tolx = ex * 1e-9, toly = ey * 1e-9, tolz = ez * 1e-9;
    std::vector<int> planes;
    if (std::fabs(p.x - cfg.min.x) < tolx) planes.push_back(0);
    if (std::fabs(p.x - cfg.max.x) < tolx) planes.push_back(1);
    if (std::fabs(p.y - cfg.min.y) < toly) planes.push_back(2);
    if (std::fabs(p.y - cfg.max.y) < toly) planes.push_back(3);
    if (std::fabs(p.z - cfg.min.z) < tolz) planes.push_back(4);
    if (std::fabs(p.z - cfg.max.z) < tolz) planes.push_back(5);
    return planes;
}

Vec3 planeNormal(int plane) {
    switch (plane) {
        case 0: return Vec3{-1, 0, 0};
        case 1: return Vec3{1, 0, 0};
        case 2: return Vec3{0, -1, 0};
        case 3: return Vec3{0, 1, 0};
        case 4: return Vec3{0, 0, -1};
        default: return Vec3{0, 0, 1};
    }
}

// --- The layer-quality GATE ------------------------------------------
// A "stack" is everything one ORIGINAL wall seed face extrudes over the
// whole march (identified by its `faceOrigin` id, minted at step 0).
// The quality evaluator above is ALWAYS on (it only touches
// freshly-created prisms, so its cost scales with the layer sheet, not
// the mesh); a prism that fails either checkMesh metric condemns its
// whole stack -- a stack that goes bad stays bad, and the damage
// concentrates at the wall-adjacent layers, so peel-back would buy
// nothing over whole-stack removal.
// Removal is expressed through the march's OWN drop path (the face is
// marked `reverted`, i.e. reverted to a plain wall face at the offset
// cut, and its seam is closed by the existing side-face rule), so the
// result is topologically identical to "this face never extruded";
// there is no post-hoc cell surgery. Because a condemned stack is only
// discovered AFTER its first prisms exist, the removal is applied by
// re-running the march with the condemned set pre-loaded, iterating to
// a fixpoint (the set only ever grows, so this terminates and is
// order-independent -- see applyLayers below).
// One condemned-stack site, opt-in dump payload.
struct CondemnedSite {
    int patchOrdinal;
    int origin;
    double x, y, z; // first-flagged prism's centroid (world coords)
    int coreLevel;  // refinement level of the core (un-extruded) cell
};

// One terraced-wall-face site, opt-in dump payload.
// A "terraced" face is any wall face that stopped marching before
// reaching its spec's full nLayers -- the SAME event the established
// `droppedFaces`/`infeasibleWallArea` stats already count (see the
// per-step drop loop in applyLayersPass), just also recorded with a
// position and level so a downstream clustering script can site the
// WHOLE degraded band, not only the (much smaller) condemned-stack
// subset. `achievedLayers` is the step index at which the face first
// dropped (0 == never completed even the first/innermost layer step);
// `specLayers` is that face's spec nLayers, so achieved-vs-specified is
// directly available without re-deriving the spec lookup downstream.
struct TerracedFaceSite {
    int patchOrdinal;
    int origin;
    double x, y, z;   // dropped face's centroid (world coords), at the offset-cut position
    int coreLevel;    // refinement level of the adjacent core (un-extruded) cell
    int achievedLayers;
    int specLayers;
    int reason;       // LayerDropReason that terraced it
};

struct LayerGateReport {
    std::set<std::pair<int, int>> badStacks; // (wall patch ordinal, origin id) flagged this pass
    // Subset of `badStacks` condemned by the SEAL guard specifically
    // (the smoothed landing crossed into the solid), kept apart so the
    // "what smoothing would have sealed" disclosure survives the
    // re-march that makes those stacks vanish from the final pass's
    // sealed-region instrument.
    std::set<std::pair<int, int>> sealedStacks;
    double sealedArea = 0.0;              // wall area of the faces the seal guard dropped this pass
    long orphanSealedFaces = 0;           // seal-guard hits on a seam face (origin -1): dropped, not condemnable
    // Same construction for the COLLAPSED-PRISM guard: stacks condemned
    // because the achieved layer height came out a small fraction of
    // the requested one. Kept apart from badStacks for the same reason
    // (the re-march makes them invisible to the final pass).
    std::set<std::pair<int, int>> thinStacks;
    double thinArea = 0.0;                // wall area of the faces the collapsed-prism guard dropped this pass
    long orphanThinFaces = 0;             // collapsed-prism hits on a seam face (origin -1)
    long totalWallFaces0 = 0;                // wall seed faces over all specs (the removal denominator)
    long orphanBadCells = 0;                 // flagged prisms grown from a seam face (origin -1): not removable
    // Opt-in (NINJA_LAYER_GATE_DUMP): one record per
    // newly-condemned stack this pass, captured at the moment its FIRST
    // bad prism is found -- the prism centroid (world coords) and the
    // core cell's refinement level, so a downstream clustering script
    // can site condemned regions without re-deriving march internals.
    // Empty (and never populated) unless the env var is set -- see
    // `gateDump` below; costs nothing on the established fast path.
    std::vector<CondemnedSite> condemnedSites;
    // Opt-in (same NINJA_LAYER_GATE_DUMP var): one
    // record per wall face that terraced (dropped short of its spec's
    // nLayers) this pass -- superset of `condemnedSites`, covering the
    // clamp-propagation collateral around the gate's own removals too.
    // Empty unless gateDump is set.
    std::vector<TerracedFaceSite> terracedFaces;
};

std::string domainPatchNameOf(int plane, const MeshConfig& cfg) {
    static const char* keys[6] = {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};
    for (const PatchSpec& ps : cfg.patches) {
        if (ps.sideKey == keys[plane]) {
            return ps.name;
        }
    }
    return "";
}

// One full march. `gateRemoved` holds the stacks condemned by earlier
// passes (never extruded here); `report` collects this pass's own
// verdicts. An empty `gateRemoved` + an empty resulting `report` is the
// pre-gate behaviour, bit for bit.
// Cumulative wall-clock per named section of one gate pass, printed at
// the end of the pass under NINJA_STAGE_TIMES (never otherwise).
struct PassClock {
    bool on = std::getenv("NINJA_STAGE_TIMES") != nullptr;
    std::chrono::steady_clock::time_point last = std::chrono::steady_clock::now();
    std::vector<std::pair<std::string, double>> acc;
    void lap(const char* name) {
        if (!on) return;
        const auto now = std::chrono::steady_clock::now();
        const double dt = std::chrono::duration<double>(now - last).count();
        last = now;
        for (auto& [n, t] : acc) {
            if (n == name) {
                t += dt;
                return;
            }
        }
        acc.emplace_back(name, dt);
    }
    void print() const {
        if (!on) return;
        double total = 0.0;
        for (const auto& [n, t] : acc) total += t;
        std::cout << "layers pass time " << total << " s:";
        for (const auto& [n, t] : acc) std::cout << " [" << n << " " << t << "]";
        std::cout << std::endl;
    }
};

LayersResult applyLayersPass(const GeneratedMesh& cutMeshIn, const std::vector<int>& cellLevelIn,
                          const std::vector<int>& pointLevelIn, const MeshConfig& cfg,
                          const std::vector<LayerStlSpec>& specs,
                          const std::vector<std::vector<Triangle>>& perStlTris,
                          const std::vector<TriangleAabbBins>& perStlBins,
                          const std::vector<ForceDropSphere>& forceDropSpheres,
                          const std::set<std::pair<int, int>>& gateRemoved,
                          LayerGateReport& report, const Vec3& locationInMesh) {
    PassClock clk; // NINJA_STAGE_TIMES: cumulative wall-clock per march section
    LayersResult result;
    GeneratedMesh out = cutMeshIn; // full copy; mutated below (points/cells appended, some
                                   // faces rewritten in place, patches rebuilt at the end)
    std::vector<int> cellLevel = cellLevelIn;
    std::vector<int> pointLevel = pointLevelIn;

    // Core cells the final flood fill will keep. The offset cut can leave
    // fluid unreachable from locationInMesh (disclosed by applyLayers); its
    // wall faces never march and never ship, so they are kept out of the
    // coverage accounting below. MEASURED, win_gate_seam: 355 of the 378
    // "dropped" faces and 89 of the 95.5 m2 of "infeasible" wall were faces
    // of one discarded 669-cell component -- none of them in the written mesh.
    // Extended with every prism as it is emitted: a prism ships iff the core
    // cell its stack grows from does, so the filter holds on EVERY step, not
    // only on the first (where the front's owner is still a core cell).
    std::vector<char> cellReached = reachableCellsFromLocation(cutMeshIn, locationInMesh);

    // name -> original patch ordinal, for locating each spec's wall patch.
    std::unordered_map<std::string, int> patchNameToOrdinal;
    for (std::size_t p = 0; p < out.patches.size(); ++p) {
        patchNameToOrdinal[out.patches[p].name] = static_cast<int>(p);
    }

    // Buckets rebuilt from scratch, in original patch order (CSR
    // FaceStores).
    FaceStore internalOut;
    for (int i = 0; i < out.nInternalFaces; ++i) {
        internalOut.appendFrom(out.faces, i);
    }
    std::vector<FaceStore> boundaryByPatch(out.patches.size());
    for (std::size_t p = 0; p < out.patches.size(); ++p) {
        const PatchInfo& pi = out.patches[p];
        for (int i = pi.startFace; i < pi.startFace + pi.nFaces; ++i) {
            boundaryByPatch[p].appendFrom(out.faces, i);
        }
    }
    // Which wall-patch ordinals (>= the domain patch count) does a layers spec cover?
    std::unordered_set<int> layeredPatchOrdinals;
    for (const LayerStlSpec& spec : specs) {
        layeredPatchOrdinals.insert(patchNameToOrdinal.at(spec.wallPatchName));
    }

    LayerStats stats;
    double residualAreaSum = 0.0;
    double residualWeightedSum = 0.0;
    double residualMax = 0.0;

    // The self-diagnosing wall-distance export needs per-node flags, not
    // just the aggregate stats above, so an artifact in
    // wallDist_<patch>.vtk can be attributed WITHOUT guessing. Grown
    // alongside `out.points` (indices are point ids in the FINAL mesh,
    // same space the point-compaction remap below operates on).
    std::vector<char> pointFrozenFlag;  // this point's origin chain was ever frozen (non-inversion clamp)
    std::vector<char> pointDroppedFlag; // this point sits on a face whose origin was ever reverted (drop-locally)

    // New prism cells are appended after the existing cell count.
    int nextCellIndex = out.nCells();
    std::vector<int> newCellLevels; // appended, in creation order
    std::vector<int> newCellStackId;     // stack topology export, parallel to newCellLevels
    std::vector<int> newCellLayerIndex;  // stack topology export, parallel to newCellLevels
    // Cell level of ANY current cell -- original cells from cellLevelIn,
    // prism cells (later layer steps stack on them) from newCellLevels.
    auto levelOfCell = [&](int c) {
        return c < static_cast<int>(cellLevelIn.size())
                   ? cellLevelIn[static_cast<std::size_t>(c)]
                   : newCellLevels[static_cast<std::size_t>(c) - cellLevelIn.size()];
    };

    // The quality evaluation itself is ALWAYS on (it is the gate's
    // predicate); `lqStats` only switches the extra REPORTING on top
    // of it.
    const bool lqStats = std::getenv("NINJA_LAYER_QUALITY_STATS") != nullptr;
    // Opt-in condemned-stack site dump, same zero-cost-unset
    // convention as NINJA_LAYER_QUALITY_STATS -- gated on its OWN var
    // so a site-dump run need not also pay for (or wade through) the
    // full quality histogram.
    const bool gateDump = std::getenv("NINJA_LAYER_GATE_DUMP") != nullptr;

    // The smoothed landing field: r = smoothRadius * h_local, per-STL
    // dict key (`LayerStlSpec::smoothRadius`). 0 disables it -- both
    // consumers below (direction, step length) then fall through to
    // their original `closestPointOnSoup` code path unconditionally, no
    // shared branch, so "byte-identical when off" does not depend on the
    // smoothed arithmetic happening to cancel out.
    // Base (level-0) cell size: r = smoothRadius * h_local, h_local =
    // dx0 / 2^level (Refine.cpp's relation). Assumes
    // cubic base cells, this project's standing convention (dx==dy==dz
    // at level 0); cheap to compute once per pass, only ever read when
    // smoothingOn.
    const double dx0 = (cfg.max.x - cfg.min.x) / std::max(1, cfg.nx);
    long lqTotalWallFaces0 = 0;                     // sum over specs of nTop at step 0
    std::set<std::pair<int, int>> lqFlaggedStacks;  // (patch ordinal, origin id)
    // Owner-side (across the top face) cell centroid lookup: for step 0
    // the owner is an ORIGINAL cutMeshIn cell (computed on demand, exact
    // polyhedron centroid, cached since several top faces can share a
    // core cell); for step > 0 the owner is the PREVIOUS step's own
    // prism, whose centroid was already computed when IT was created.
    std::unordered_map<int, Vec3> lqOrigCellCentroidCache;
    std::unordered_map<int, Vec3> lqPrismCentroid; // prism cell id -> its own centroid
    auto lqOwnerCentroidOf = [&](int coreCell) -> const Vec3* {
        if (coreCell < cutMeshIn.nCells()) {
            auto it = lqOrigCellCentroidCache.find(coreCell);
            if (it != lqOrigCellCentroidCache.end()) return &it->second;
            const Vec3 c = cellCentroidOf(cutMeshIn, coreCell);
            return &lqOrigCellCentroidCache.emplace(coreCell, c).first->second;
        }
        auto it = lqPrismCentroid.find(coreCell);
        return it != lqPrismCentroid.end() ? &it->second : nullptr;
    };

    clk.lap("pass setup");
    for (const LayerStlSpec& spec : specs) {
        const int pOrd = patchNameToOrdinal.at(spec.wallPatchName);
        FaceStore& wallBucket = boundaryByPatch[static_cast<std::size_t>(pOrd)];
        const TriangleAabbBins& bins = perStlBins[static_cast<std::size_t>(spec.stlIndex)];
        const double t = spec.thickness;
        // Per-spec smoothing radius fraction -- env wins over the
        // dict when both are present (documented precedence, `layers{}`
        // dict comment).
        const double smoothRadius = spec.smoothRadius;
        const bool smoothingOn = smoothRadius > 0.0;
        const SmoothedField smoothCfg{smoothRadius};

        // --- Layer schedule (layer-by-layer march):
        // nLayers graded steps, thinnest layer AT THE WALL, marched
        // outermost-first. Layer j (j = 0 at the wall) has thickness
        // t0 * ratio^j with sum == t. Step lengths are expressed in
        // FIELD units: at each step the target advance is that layer's
        // fraction of the point's own REMAINING distance d(p), so
        // points clamped short on an earlier step automatically re-aim
        // with what distance they have left, and the last step (the
        // wall layer) always targets the full remaining d.
        const int nSteps = std::max(1, spec.nLayers);
        const double ratio = spec.expansionRatio > 0.0 ? spec.expansionRatio : 1.0;
        std::vector<double> layerThickness(static_cast<std::size_t>(nSteps));
        {
            double denom = 0.0;
            for (int j = 0; j < nSteps; ++j) denom += std::pow(ratio, j);
            const double t0 = t / denom;
            for (int j = 0; j < nSteps; ++j) layerThickness[static_cast<std::size_t>(j)] = t0 * std::pow(ratio, j);
        }

        // --- DISTINCT drop/freeze bookkeeping.
        // `faceOrigin` runs strictly parallel to `wallBucket` (never
        // reordered; rebuilt in the same push order as newWallBucket)
        // and carries, for each current wall face, the index of the
        // ORIGINAL layer-top face it descends from -- or -1 for faces
        // that have no such ancestor (seam side faces created by a
        // drop). `pointOrigin` does the same for front points. A face
        // (point) is added to `stats` only the first time its origin
        // enters the dropped (frozen) set, so a face reverted on
        // several successive steps counts once.
        std::vector<int> faceOrigin(static_cast<std::size_t>(wallBucket.size()));
        for (std::size_t i = 0; i < faceOrigin.size(); ++i) faceOrigin[i] = static_cast<int>(i);
        // A front point is extruded AT MOST ONCE. Once a kept prism has
        // marched point P to P', P is spent: a face still holding P (a
        // terrace seam, or a dropped face beside the kept one) is HELD --
        // it stays a wall face for the rest of the march. Were it to
        // extrude later, P would march a second time down the same line
        // as P' and both land on the same wall spot: duplicate points and
        // a zero-width crack of back-to-back wall faces (MEASURED,
        // bm_layers_gate: 9 duplicate positions). Seam faces satisfy this
        // by construction (their top edge is spent).
        std::vector<char> pointSpent;
        auto faceHeld = [&](int fi) {
            if (faceOrigin[static_cast<std::size_t>(fi)] < 0) return true;
            for (int p : wallBucket.pointsOf(fi)) {
                if (static_cast<std::size_t>(p) < pointSpent.size() && pointSpent[static_cast<std::size_t>(p)]) return true;
            }
            return false;
        };
        // The STACK identity the quality gate condemns
        // by. Starts out equal to `faceOrigin` and is inherited the same
        // way, with ONE deliberate difference: a terrace SEAM face
        // (created by the kept prism `fi` alongside a dropped/removed
        // neighbour, `faceOrigin` -1 because it has no wall ancestor)
        // inherits the identity of the prism that CREATED it. MEASURED
        // necessary (offset_rotcube): seam faces are re-extruded on the
        // next step, those prisms can themselves be degenerate, and with
        // no stack to condemn the gate could not converge -- it left 16
        // bad cells (8 collinear, zero-area side quads -> 8 open cells)
        // in the mesh. Attributing a bad seam prism to the stack whose
        // extrusion produced the seam removes that stack, so the seam
        // never exists. Kept SEPARATE from `faceOrigin` so the
        // drop/freeze accounting (droppedFaces, infeasibleWallArea,
        // pointDropped) keeps its own established meaning.
        std::vector<int> faceGateStack = faceOrigin;
        std::unordered_map<int, int> pointOrigin; // point index -> origin id
        std::set<int> droppedOrigins;
        // Every layer prism by cell: the loops it was built between and the
        // cell it grew from -- what a prism merged into it (see mergeInto)
        // needs. A prism absorbs at most one failed prism below it.
        std::unordered_map<int, std::vector<int>> prismTopOf, prismBottomOf;
        std::unordered_map<int, int> prismParentOf;
        std::unordered_set<int> mergedPrisms;
        std::set<int> frozenOrigins;
        // Quality instrumentation: origin id -> wall seed-face non-planarity,
        // captured once at step 0 (where wallBucket IS the original
        // cutter-emitted wall geometry, before any marching) since
        // origin ids are minted there (faceOrigin[fi] = fi, above).
        std::unordered_map<int, double> lqSeedNonPlanarity;
        auto originOfPoint = [&](int p) {
            auto it = pointOrigin.find(p);
            if (it != pointOrigin.end()) return it->second;
            // First sighting: the point is its own origin (point indices
            // are globally unique, so this can never alias).
            pointOrigin.emplace(p, p);
            return p;
        };
        if (static_cast<int>(stats.perStepDropped.size()) < nSteps) {
            stats.perStepDropped.resize(static_cast<std::size_t>(nSteps), 0);
            stats.perStepFrozen.resize(static_cast<std::size_t>(nSteps), 0);
            stats.perStepPrismCells.resize(static_cast<std::size_t>(nSteps), 0);
        }

        clk.lap("spec setup");
        for (int step = 0; step < nSteps; ++step) {
        const int layerJ = nSteps - 1 - step; // marching outermost layer first
        const double tStep = layerThickness[static_cast<std::size_t>(layerJ)];
        double remThickness = 0.0; // this layer + everything still inside it
        for (int j = 0; j <= layerJ; ++j) remThickness += layerThickness[static_cast<std::size_t>(j)];
        const bool lastStep = layerJ == 0;
        // Field-unit fraction of the remaining distance this step targets
        // (exactly 1.0 on the wall layer).
        const double stepFrac = lastStep ? 1.0 : tStep / remThickness;
        const double minVolEps = 1e-9 * tStep * tStep * tStep;
        const double minHeight = kMinHeightFrac * tStep;

        const int nTop = wallBucket.size();
        if (step == 0) lqTotalWallFaces0 += nTop;
        // Quality-instrumentation per-step accumulators (cheap: only populated
        // when lqStats is set; default-constructed containers below
        // never allocate otherwise).
        int lqBadPyramidFaces = 0, lqBadTetFaces = 0, lqFlaggedCells = 0;
        std::vector<double> lqWorstPerCell;
        std::vector<std::pair<int, double>> lqNewFlaggedThisStep; // (origin id, non-planarity)

        clk.lap("step header");
        // --- Front construction -------------------------------------
        std::set<int> frontSet;
        for (int fi = 0; fi < nTop; ++fi) {
            for (int p : wallBucket.pointsOf(fi)) frontSet.insert(p);
        }
        std::vector<int> frontPoints(frontSet.begin(), frontSet.end());
        std::unordered_map<int, int> frontIdx; // orig point idx -> dense index
        for (std::size_t i = 0; i < frontPoints.size(); ++i) frontIdx[frontPoints[i]] = static_cast<int>(i);
        const std::size_t nFront = frontPoints.size();

        std::vector<std::unordered_set<int>> adjacency(nFront); // dense-index adjacency
        std::unordered_map<EdgeKey, std::vector<int>, EdgeKeyHash> edgeToFaces; // face local index
        for (int fi = 0; fi < nTop; ++fi) {
            const IntSpan fp = wallBucket.pointsOf(fi);
            const int n = fp.size();
            for (int i = 0; i < n; ++i) {
                const int a = fp[i];
                const int b = fp[(i + 1) % n];
                adjacency[static_cast<std::size_t>(frontIdx[a])].insert(frontIdx[b]);
                adjacency[static_cast<std::size_t>(frontIdx[b])].insert(frontIdx[a]);
                edgeToFaces[makeEdgeKey(a, b)].push_back(fi);
            }
        }

        // `localThickness`: each face's thickness is the local field at the
        // face -- faceScale = t(face) / t (1 otherwise); a point takes the
        // thinnest of the faces around it. Only the ABSOLUTE thresholds
        // scale: the step fractions below are ratios of layer
        // thicknesses, and the front already sits at the local offset
        // distance the cut gave it (from the same field).
        std::vector<double> faceScale(static_cast<std::size_t>(nTop), 1.0);
        std::vector<double> pointScale(nFront, 1.0);
        if (spec.localThickness) {
            std::fill(pointScale.begin(), pointScale.end(), std::numeric_limits<double>::max());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256)
#endif
            for (int fi = 0; fi < nTop; ++fi) {
                std::vector<Vec3> loop;
                for (int p : wallBucket.pointsOf(fi)) loop.push_back(out.points[static_cast<std::size_t>(p)]);
                faceScale[static_cast<std::size_t>(fi)] = spec.localThickness(loopCentroid(loop)) / t;
            }
            for (int fi = 0; fi < nTop; ++fi) {
                const double sc = faceScale[static_cast<std::size_t>(fi)];
                for (int p : wallBucket.pointsOf(fi)) {
                    double& ps = pointScale[static_cast<std::size_t>(frontIdx[p])];
                    ps = std::min(ps, sc);
                }
            }
        }

        std::vector<double> dVal(nFront, 0.0);
        std::vector<Vec3> hitPoint(nFront);
        std::vector<Vec3> dir(nFront);
        std::vector<Vec3> frontNormal(nFront, Vec3{0, 0, 0});
        std::vector<std::vector<int>> domainPlanes(nFront);
        // Per-point smoothing radius, r = smoothRadius * h_local --
        // computed (and kept, for the landing-time sealed-region
        // instrumentation below) regardless of smoothingOn's value would
        // be wasted work when off, so this is itself gated.
        std::vector<double> radii;
        if (smoothingOn) {
            radii.assign(nFront, 0.0);
            for (int fi = 0; fi < nTop; ++fi) {
                const int lvl = levelOfCell(wallBucket.owner[static_cast<std::size_t>(fi)]);
                const double h = dx0 / static_cast<double>(1 << lvl);
                for (int p : wallBucket.pointsOf(fi)) {
                    const std::size_t idx = static_cast<std::size_t>(frontIdx[p]);
                    // The FINEST level touching the point decides: at a
                    // level transition the coarse side's radius would smooth
                    // the fine side's wall over twice its own cell size.
                    radii[idx] = radii[idx] > 0.0 ? std::min(radii[idx], smoothRadius * h) : smoothRadius * h;
                }
            }
        }

        if (!smoothingOn) {
            // --- Exact-field path. The flag being unset takes this
            // path unconditionally, which is the byte-identity
            // argument: smoothing off never touches smoothed
            // arithmetic at all.
            // Pure per-point aim; parallel (each iteration writes only its
            // own slots of dVal/hitPoint/dir/domainPlanes).
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 512)
#endif
            for (std::ptrdiff_t fi = 0; fi < static_cast<std::ptrdiff_t>(nFront); ++fi) {
                const std::size_t i = static_cast<std::size_t>(fi);
                const Vec3& pos = out.points[static_cast<std::size_t>(frontPoints[i])];
                ClosestHit hit = closestPointOnSoup(bins, pos);
                dVal[i] = std::sqrt(hit.distSq);
                hitPoint[i] = hit.point;
                // March direction: DESCENT along the distance field, i.e.
                // toward the closest surface point (-grad d), NOT
                // normalize(p - closestPointOnSoup(p)) (which is +grad d,
                // AWAY from the surface) -- the method only makes sense
                // stepping the front DOWN toward d=0. Treated as a corrected sign,
                // implementing -grad d = normalize(closest - p).
                dir[i] = dVal[i] > 1e-300 ? normalizeOrZero(hit.point - pos) : Vec3{0, 0, 0};
                domainPlanes[i] = domainPlanesOf(pos, cfg);
            }
        } else {
            // --- Smoothed path: BOTH march consumers -- direction and
            // the remaining-distance used for step length -- read off
            // the SMOOTHED field, batched over the whole front in one
            // call (one classifyVertices parity pass for every stencil
            // point this step needs).
            std::vector<Vec3> qpts(nFront);
            for (std::size_t i = 0; i < nFront; ++i) qpts[i] = out.points[static_cast<std::size_t>(frontPoints[i])];
            std::vector<Vec3> gradPhi;
            const std::vector<double> phi = smoothedSignedDistanceBatch(
                bins, perStlTris[static_cast<std::size_t>(spec.stlIndex)], locationInMesh, qpts, radii, smoothCfg,
                &gradPhi);
            for (std::size_t i = 0; i < nFront; ++i) {
                dVal[i] = phi[i]; // signed remaining distance in field units (may be <= 0 near a sealed feature)
                hitPoint[i] = qpts[i]; // no exact closest point in this path; kept only for shape/debug parity
                // Direction: normalize(-grad phi_r), read directly off
                // the same stencil already evaluated for phi_r itself
                // (Geometry.hpp's smoothedSignedDistanceBatch gradOut,
                // no extra queries) -- NOT closestPoint - p (that would
                // silently fall back to the exact field's discontinuity,
                // exactly what smoothing exists to remove). A
                // degenerate estimate (|grad| ~ 0) is counted, never
                // silently substituted, and falls through to the
                // existing frontNormal fallback below (the SAME fallback
                // the exact path already uses for its own degenerate
                // case).
                const double gnorm = norm(gradPhi[i]);
                if (gnorm > 1e-300) {
                    dir[i] = gradPhi[i] * (-1.0 / gnorm);
                } else {
                    dir[i] = Vec3{0, 0, 0};
                    ++stats.smoothGradFallback;
                }
                domainPlanes[i] = domainPlanesOf(qpts[i], cfg);
            }
        }
        // Area-weighted front normal (outward-from-core, as stored).
        for (int fi = 0; fi < nTop; ++fi) {
            const IntSpan fp = wallBucket.pointsOf(fi);
            std::vector<Vec3> loop;
            loop.reserve(static_cast<std::size_t>(fp.size()));
            for (int p : fp) loop.push_back(out.points[static_cast<std::size_t>(p)]);
            const Vec3 nw = newellNormal(loop); // magnitude ~ 2*area, direction outward
            for (int p : fp) {
                frontNormal[static_cast<std::size_t>(frontIdx[p])] = frontNormal[static_cast<std::size_t>(frontIdx[p])] + nw;
            }
        }
        for (std::size_t i = 0; i < nFront; ++i) frontNormal[i] = normalizeOrZero(frontNormal[i]);
        // Points with a degenerate initial direction (sitting exactly on
        // their closest point) march along -frontNormal instead.
        for (std::size_t i = 0; i < nFront; ++i) {
            if (norm(dir[i]) < 1e-300) {
                dir[i] = frontNormal[i]; // descent direction fallback, see refDir note below
            }
        }

        clk.lap("front + field");
        // --- Laplacian smoothing, weights favouring similar remaining d
        // (documented formula: w(p,q) = 1 / (1 + |d(p)-d(q)| / t)).
        for (int pass = 0; pass < kLaplacianPasses; ++pass) {
            std::vector<Vec3> next = dir;
            for (std::size_t i = 0; i < nFront; ++i) {
                if (adjacency[i].empty()) continue;
                Vec3 acc{0, 0, 0};
                double wsum = 0.0;
                for (int j : adjacency[i]) {
                    const double w = 1.0 / (1.0 + std::fabs(dVal[i] - dVal[static_cast<std::size_t>(j)]) / std::max(t * pointScale[i], 1e-300));
                    acc = acc + dir[static_cast<std::size_t>(j)] * w;
                    wsum += w;
                }
                if (wsum <= 0.0) continue;
                const Vec3 nbrAvg = acc * (1.0 / wsum);
                const Vec3 blended = dir[i] * kLaplacianBlend + nbrAvg * (1.0 - kLaplacianBlend);
                next[i] = normalizeOrZero(blended);
                if (norm(next[i]) < 1e-300) next[i] = dir[i];
            }
            dir = std::move(next);
        }

        clk.lap("laplacian");
        // --- Domain-boundary in-plane projection.
        for (std::size_t i = 0; i < nFront; ++i) {
            if (domainPlanes[i].empty()) continue;
            Vec3 d = dir[i];
            for (int plane : domainPlanes[i]) {
                const Vec3 n = planeNormal(plane);
                d = d - n * dot(d, n);
            }
            dir[i] = normalizeOrZero(d); // zero at a fully-constrained corner: point cannot move in-plane
        }

        clk.lap("domain projection");
        // --- Micro-edge rigid clusters (see kMicroEdgeFrac): union-find
        // over front edges shorter than kMicroEdgeFrac * tStep, then
        // unify the march direction across each cluster (normalized
        // average; step lengths are equalized after the clamps below).
        std::vector<int> clusterOf(nFront);
        for (std::size_t i = 0; i < nFront; ++i) clusterOf[i] = static_cast<int>(i);
        {
            std::function<int(int)> findC = [&](int x) {
                while (clusterOf[static_cast<std::size_t>(x)] != x) {
                    clusterOf[static_cast<std::size_t>(x)] =
                        clusterOf[static_cast<std::size_t>(clusterOf[static_cast<std::size_t>(x)])];
                    x = clusterOf[static_cast<std::size_t>(x)];
                }
                return x;
            };
            const double microLen = kMicroEdgeFrac * tStep;
            for (std::size_t i = 0; i < nFront; ++i) {
                const Vec3& p0 = out.points[static_cast<std::size_t>(frontPoints[i])];
                for (int j : adjacency[i]) {
                    if (static_cast<std::size_t>(j) <= i) continue;
                    const Vec3& q0 = out.points[static_cast<std::size_t>(frontPoints[static_cast<std::size_t>(j)])];
                    if (norm(p0 - q0) < microLen * std::min(pointScale[i], pointScale[static_cast<std::size_t>(j)])) {
                        const int ra = findC(static_cast<int>(i));
                        const int rb = findC(j);
                        if (ra != rb) clusterOf[static_cast<std::size_t>(ra)] = rb;
                    }
                }
            }
            for (std::size_t i = 0; i < nFront; ++i) clusterOf[i] = findC(static_cast<int>(i));
            // Unify dir (and the angle-clamp reference normal) per cluster.
            std::vector<Vec3> dirAcc(nFront, Vec3{0, 0, 0}), nrmAcc(nFront, Vec3{0, 0, 0});
            std::vector<int> cnt(nFront, 0);
            for (std::size_t i = 0; i < nFront; ++i) {
                const std::size_t r = static_cast<std::size_t>(clusterOf[i]);
                dirAcc[r] = dirAcc[r] + dir[i];
                nrmAcc[r] = nrmAcc[r] + frontNormal[i];
                ++cnt[r];
            }
            for (std::size_t i = 0; i < nFront; ++i) {
                const std::size_t r = static_cast<std::size_t>(clusterOf[i]);
                if (cnt[r] > 1) {
                    dir[i] = normalizeOrZero(dirAcc[r]);
                    frontNormal[i] = normalizeOrZero(nrmAcc[r]);
                }
            }
        }

        clk.lap("micro-edge clusters");
        // --- Angle clamp -> per-point UNCLAMPED-for-stretch step length.
        // frontNormal is the area-weighted average of the STORED
        // top-face orientation, which is outward-FROM-CORE -- and since
        // the core (fluid) cell sits on the d > t side of the offset
        // boundary, "outward from core" already points in the DESCENT
        // direction (toward smaller d, toward the true STL). So the
        // march direction should align WITH frontNormal, not against
        // it.
        std::vector<double> stepLen0(nFront, 0.0);
        for (std::size_t i = 0; i < nFront; ++i) {
            if (norm(dir[i]) < 1e-300 || dVal[i] <= 0.0) {
                continue;
            }
            double scaleAngle = 1.0;
            const Vec3 refDir = frontNormal[i];
            if (norm(refDir) > 1e-300) {
                const double cosA = std::clamp(dot(dir[i], refDir), -1.0, 1.0);
                const double angleDeg = std::acos(cosA) * 180.0 / M_PI;
                if (angleDeg > kMaxAngleDeg) {
                    scaleAngle = std::clamp(kMaxAngleDeg / angleDeg, 0.0, 1.0);
                }
            }
            stepLen0[i] = dVal[i] * stepFrac * scaleAngle;
        }

        // Tangential stretch: a whole-front front-graph is converging
        // (or expanding) roughly IN SYNC (every point marches on the
        // same step), so the representative candidate for "how does
        // this edge change" is BOTH endpoints advancing together by the
        // same fraction of their own (already angle-clamped) step --
        // not one endpoint moving while its neighbour stays pinned at
        // its ORIGINAL (offset) position, which manufactures spurious
        // stretch/compression for any front that is simply
        // contracting uniformly toward a smaller true surface (e.g. a
        // sphere: EVERY edge shortens by the same modest ratio, well
        // inside kStretchMax, but measuring one moving endpoint against
        // a static neighbour makes it look enormous). Found the largest
        // common fraction `frac` (of each point's OWN stepLen0) via a
        // monotonic-safe forward scan (K samples, per-edge), local and
        // final -- no global relax-and-retry loop.
        std::vector<double> frac(nFront, 1.0);
        constexpr int kStretchSamples = 32;
        for (std::size_t i = 0; i < nFront; ++i) {
            if (stepLen0[i] <= 0.0) { frac[i] = 0.0; continue; }
            const Vec3& p0 = out.points[static_cast<std::size_t>(frontPoints[i])];
            for (int j : adjacency[i]) {
                const std::size_t jz = static_cast<std::size_t>(j);
                // Intra-cluster edges march rigidly -- no relative
                // stretch by construction, and their micro origLen would
                // otherwise blow the ratio up and crush frac to ~0 (the
                // measured seed of the sliver-facet failures).
                if (clusterOf[i] == clusterOf[jz]) continue;
                const Vec3& q0 = out.points[static_cast<std::size_t>(frontPoints[jz])];
                const double origLen = norm(p0 - q0);
                if (origLen < 1e-300) continue;
                auto ratioAt = [&](double f) {
                    const Vec3 pCand = p0 + dir[i] * (stepLen0[i] * f);
                    const Vec3 qCand = q0 + dir[jz] * (stepLen0[jz] * f);
                    return norm(pCand - qCand) / origLen;
                };
                double lastGood = 0.0;
                for (int k = 1; k <= kStretchSamples; ++k) {
                    const double f = static_cast<double>(k) / kStretchSamples;
                    const double r = ratioAt(f);
                    if (r <= kStretchMax && r >= 1.0 / kStretchMax) {
                        lastGood = f;
                    } else {
                        break; // first violation walking forward from f=0 (always valid, r=1)
                    }
                }
                frac[i] = std::min(frac[i], lastGood);
            }
        }
        std::vector<double> stepLen(nFront, 0.0);
        for (std::size_t i = 0; i < nFront; ++i) {
            stepLen[i] = stepLen0[i] * frac[i];
        }
        clk.lap("angle clamp");
        // --- Surface guard: a step may not carry a front point THROUGH a
        // surface and leave it in the fluid on the far side. A sound landing can
        // pass the true surface -- the smoothed field rounds convex edges, and
        // its zero lies just inside the solid there -- but then parity says
        // SOLID at the new position. A point that crosses the tangent plane at
        // its exact closest point and is still FLUID by parity has gone through
        // a wall the signed field cannot see; it is landed on that wall instead.
        // MEASURED, win_desc_slot: a second CAD body's face lies exactly on the
        // pier's outer face with the opposite normal. The offset cut measures
        // unsigned distance and builds the front there; in the signed field two
        // coincident crossings cancel, so it is a V with no sign change whose
        // smoothed bottom never reaches zero. The front stepped through the
        // wall (third interface 14 mm past it), the landing step turned round,
        // and the fold test refused 1414 prisms -- by luck: clamping only the
        // intermediate steps let the landing run 35 mm into the pier as a
        // perfectly valid prism. On sound CAD it binds a handful of times per
        // case, at plates thinner than a step (1 on win_fish_coarse_site1, 10 on
        // win_desc_coarse), and 7 of 9 layered windows come out byte-identical.
        if (smoothingOn) {
            std::vector<ClosestHit> exactHit(nFront);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 512)
#endif
            for (std::ptrdiff_t fi = 0; fi < static_cast<std::ptrdiff_t>(nFront); ++fi) {
                const std::size_t i = static_cast<std::size_t>(fi);
                if (stepLen[i] > 0.0) exactHit[i] = closestPointOnSoup(bins, out.points[static_cast<std::size_t>(frontPoints[i])]);
            }
            std::vector<std::size_t> crossing;
            std::vector<Vec3> crossedTo;
            for (std::size_t i = 0; i < nFront; ++i) {
                if (stepLen[i] <= 0.0 || stepLen[i] * stepLen[i] <= exactHit[i].distSq) continue;
                const Vec3& p0 = out.points[static_cast<std::size_t>(frontPoints[i])];
                const Vec3 pNew = p0 + dir[i] * stepLen[i];
                // Beyond the tangent plane at the closest point (a point already
                // ON the surface has no side left to stay on).
                if (dot(pNew - exactHit[i].point, p0 - exactHit[i].point) > 0.0) continue;
                crossing.push_back(i);
                crossedTo.push_back(pNew);
            }
            if (!crossing.empty()) {
                const std::vector<bool> solidThere =
                    classifyVertices(crossedTo, perStlTris[static_cast<std::size_t>(spec.stlIndex)], locationInMesh, &bins);
                for (std::size_t k = 0; k < crossing.size(); ++k) {
                    if (solidThere[k]) continue;
                    stepLen[crossing[k]] = std::sqrt(exactHit[crossing[k]].distSq);
                    ++stats.surfaceClampedSteps;
                }
            }
        }

        // Rigid clusters share ONE step length (the cluster minimum, so
        // every member honours every member's clamps).
        {
            std::vector<double> clusterMin(nFront, std::numeric_limits<double>::max());
            for (std::size_t i = 0; i < nFront; ++i) {
                const std::size_t r = static_cast<std::size_t>(clusterOf[i]);
                clusterMin[r] = std::min(clusterMin[r], stepLen[i]);
            }
            for (std::size_t i = 0; i < nFront; ++i) {
                stepLen[i] = clusterMin[static_cast<std::size_t>(clusterOf[i])];
            }
        }

        std::vector<Vec3> tentative(nFront);
        for (std::size_t i = 0; i < nFront; ++i) {
            tentative[i] = out.points[static_cast<std::size_t>(frontPoints[i])] + dir[i] * stepLen[i];
        }
        clk.lap("surface guard");
        // --- Non-inversion clamp: bisect ALL points of any invalid
        // prism toward their original (top) position, iterate, freeze
        // if still invalid after the bound.
        //
        // CASCADE FIX (measured on offset_sphere + merged slivers,
        // VTK-instrumented): the bisection must NOT enforce
        // the minHeight floor (prismValid is called with minHeight 0
        // here). Halving a point SHRINKS its prism heights, so a face
        // whose only defect is "an already-halved neighbour point made
        // an edge shorter than minHeight" would fail too, halving ITS
        // other points -- freezing then spreads one ring per outer
        // iteration and eats the whole front (measured: 63 seed
        // failures -> 150/150 points frozen, all 152 faces dropped).
        // Bisection now reacts only to genuine twist/inversion (which
        // halving actually improves); the minHeight floor is enforced
        // at FINAL validation below, where a too-thin prism is dropped
        // LOCALLY instead of contaminating its neighbours.
        std::vector<double> scaleFactor(nFront, 1.0);
        auto botPosOf = [&](std::size_t i) {
            const Vec3& top = out.points[static_cast<std::size_t>(frontPoints[i])];
            return top + (tentative[i] - top) * scaleFactor[i];
        };
        // Convergence guard on the non-inversion clamp (MEASURED, real
        // bathymetry): the bisection does NOT converge there -- the
        // invalid-face count plateaus after ~2 iterations (step 5 of a
        // 80x80 m bathymetry region: 229 invalid at outer 0, then
        // 221/212/208/216/213/212/212/210/212) while every point it
        // touches keeps halving to the 2^-10 cap. Worse, halving is
        // actively counterproductive for the "bottom face folded"
        // failure (code 6: 42 -> 82 over those same ten iterations) --
        // collapsing the bottom loop onto the top loop is exactly what
        // makes its normal degenerate.
        //
        // The crushed points then miss the minHeight floor at FINAL
        // validation, so their faces are dropped and terraced, which
        // enlarges the front and feeds the next step: a cascade from
        // 459 marched prisms at step 1 down to 87 at step 5, and 90% of
        // the wall area reported infeasible.
        //
        // So: keep the halving that demonstrably helps and stop the
        // moment it stops helping, restoring the best scale seen. This
        // introduces no tolerance -- the criterion is the clamp's own
        // invalid-face count -- and leaves what is still invalid to be
        // dropped LOCALLY by final validation, which is what the
        // cascade-fix comment above already intends. Measured on that
        // region: infeasible wall area 89.7% -> 44.6%, per-step prism
        // cells 408/256/117/41/8 -> 472/422/353/296/238, frozen points
        // -> 0.
        std::vector<double> bestScale = scaleFactor;
        int bestInvalid = std::numeric_limits<int>::max();
        for (int outer = 0; outer < kNonInversionMaxIter; ++outer) {
            bool anyInvalid = false;
            int nInvalidThisIter = 0;
            // At most ONE halving per point per outer iteration -- a
            // point touching several invalid faces in the same sweep
            // must not get compounded reductions (that turns one
            // genuinely-thin face into a runaway cascade across the
            // whole front, since a lopsided-scale corner can then make
            // an otherwise-fine neighbouring face look invalid too).
            std::vector<bool> haltThisIter(nFront, false);
            for (int fi = 0; fi < nTop; ++fi) {
                const IntSpan fp = wallBucket.pointsOf(fi);
                std::vector<Vec3> topLoop, botLoop;
                topLoop.reserve(static_cast<std::size_t>(fp.size()));
                botLoop.reserve(static_cast<std::size_t>(fp.size()));
                for (int p : fp) {
                    topLoop.push_back(out.points[static_cast<std::size_t>(p)]);
                    botLoop.push_back(botPosOf(static_cast<std::size_t>(frontIdx[p])));
                }
                // A held face never extrudes, so its would-be prism must not
                // drive the clamp: a seam fails by construction (the face is
                // perpendicular to the march) and halving its points crushes
                // the healthy stacks that share them.
                if (faceHeld(fi)) continue;
                const double fs = faceScale[static_cast<std::size_t>(fi)];
                const bool valid = prismValid(topLoop, botLoop, /*minHeight=*/0.0, minVolEps * fs * fs * fs);
                if (!valid) {
                    anyInvalid = true;
                    ++nInvalidThisIter;
                    for (int p : fp) haltThisIter[static_cast<std::size_t>(frontIdx[p])] = true;
                }
            }
            if (!anyInvalid) break;
            if (nInvalidThisIter < bestInvalid) {
                bestInvalid = nInvalidThisIter;
                bestScale = scaleFactor;
            } else {
                // This round of halving did not shrink the invalid set:
                // going further only crushes the front. Restore the best
                // scale seen and let final validation drop, locally,
                // whatever is still invalid.
                scaleFactor = bestScale;
                break;
            }
            // Halving is cluster-synchronized: if any member of a rigid
            // cluster halves, the whole cluster halves (otherwise the
            // rigidity that removes micro-edge distortion is lost the
            // moment the non-inversion clamp touches one member).
            for (std::size_t i = 0; i < nFront; ++i) {
                if (haltThisIter[i]) {
                    haltThisIter[static_cast<std::size_t>(clusterOf[i])] = true;
                }
            }
            for (std::size_t i = 0; i < nFront; ++i) {
                if (haltThisIter[static_cast<std::size_t>(clusterOf[i])]) scaleFactor[i] *= 0.5;
            }
        }
        for (std::size_t i = 0; i < nFront; ++i) {
            if (scaleFactor[i] < kFreezeScaleThreshold) {
                ++stats.perStepFrozen[static_cast<std::size_t>(step)];
                if (frozenOrigins.insert(originOfPoint(frontPoints[i])).second) ++stats.frozenPoints;
            }
        }
        clk.lap("non-inversion clamp");
        // --- Landing (LAST layer step only): re-query at the (possibly
        // frozen/clamped) tentative position. Intermediate steps must
        // NOT snap/land -- their fronts are meant to sit at the layer
        // interfaces, not the wall; the snap/landing threshold is sized
        // against the WALL layer's own thickness either way.
        //
        // smoothingOff: snap to the exact closest point on the true STL
        // if nearly there.
        //
        // smoothingOn: land on the phi_r=0 level set instead. Landing on
        // the exact surface while marching on the smoothed field would
        // mix the two fields and reintroduce the landing error the
        // smoothed path exists to avoid. The root find must be CHEAP:
        // every field evaluation pays a full batched
        // smoothedSignedDistanceBatch call (29 closest-point queries +
        // one classifyVertices parity pass per point), and a naive
        // bracket-doubling + bisection scheme spends ~80 evaluations per
        // landing candidate over hundreds of thousands of stacks --
        // enough to dominate the whole mesher. The solver below exploits
        // that phi_r is very nearly linear along the descent direction
        // (|grad phi_r| ~ 1 away from creases), so it converges in a few
        // evaluations, with a bounded bracketed fallback. Every fallback
        // (down to the exact snap) is counted, never silent -- a silent
        // substitution here would reintroduce the landing error
        // invisibly.
        std::vector<Vec3> finalPos(nFront);
        std::vector<Vec3> preLandAll(nFront);
        for (std::size_t i = 0; i < nFront; ++i) preLandAll[i] = botPosOf(i);

        std::vector<char> nearlyThere(nFront, 0);
        std::vector<Vec3> exactSnap(nFront);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 512)
#endif
        for (std::ptrdiff_t fi = 0; fi < static_cast<std::ptrdiff_t>(nFront); ++fi) {
            const std::size_t i = static_cast<std::size_t>(fi);
            if (lastStep && scaleFactor[i] >= kFreezeScaleThreshold) {
                ClosestHit hit2 = closestPointOnSoup(bins, preLandAll[i]);
                const double remaining = std::sqrt(hit2.distSq);
                nearlyThere[i] = remaining < kLandFrac * tStep * pointScale[i] ? 1 : 0;
                exactSnap[i] = hit2.point;
            }
        }

        for (std::size_t i = 0; i < nFront; ++i) {
            if (!(lastStep && scaleFactor[i] >= kFreezeScaleThreshold) || !nearlyThere[i]) {
                finalPos[i] = preLandAll[i]; // frozen, intermediate, or not yet near: no snap/land
            } else if (!smoothingOn) {
                finalPos[i] = exactSnap[i];
            }
            // else: smoothingOn landing candidates are resolved below,
            // batched over the whole candidate set.
        }

        if (smoothingOn) {
            std::vector<std::size_t> candIdx;
            for (std::size_t i = 0; i < nFront; ++i) {
                if (lastStep && scaleFactor[i] >= kFreezeScaleThreshold && nearlyThere[i]) candIdx.push_back(i);
            }
            if (!candIdx.empty()) {
                const std::size_t nc = candIdx.size();
                std::vector<Vec3> pts0(nc);
                std::vector<double> r0(nc);
                for (std::size_t k = 0; k < nc; ++k) {
                    pts0[k] = preLandAll[candIdx[k]];
                    r0[k] = radii[candIdx[k]];
                }
                std::vector<Vec3> grad0(nc);
                const std::vector<double> phi0 = smoothedSignedDistanceBatch(
                    bins, perStlTris[static_cast<std::size_t>(spec.stlIndex)], locationInMesh, pts0, r0, smoothCfg,
                    &grad0);

                // Solver choice: Illinois-safeguarded regula falsi -- a
                // BRACKETING method (can never step outside [lo, hi] and
                // hence needs no divergence-rescue branch or failure
                // counter) that, because phi_r is very nearly LINEAR
                // along the descent direction (|grad phi_r| ~ 1 by
                // construction of a smoothed distance field away from
                // creases), converges about as fast as Newton in
                // practice. Newton/secant (reusing the batched gradient)
                // is deliberately NOT used: it is not a bracketing
                // method, so a bad step near a crease/medial-axis kink
                // can overshoot arbitrarily and silently change the
                // ANSWER (drop/feasibility results), not just the cost.
                // Plain false position is also avoided: it stalls when
                // one bracket endpoint is repeatedly retained (linear,
                // not superlinear, convergence); the Illinois safeguard
                // (halve the STICKING endpoint's stored function value
                // whenever it is retained twice in a row) restores
                // superlinear convergence with a few lines and is kept
                // even though it looks removable -- deleting it
                // reintroduces the stall.
                //
                // Bracket seeding: phi0 must NOT be assumed positive
                // (fluid) even though `nearlyThere` gates on the EXACT
                // closest-point distance -- the SMOOTHED field can
                // already be negative at the pre-landing sample (the
                // rounding has "eaten" the local wall by more than the
                // exact remaining distance; expected, this is the
                // achieved-thickness confound). Marching along -grad
                // phi_r from a point already on the solid side moves
                // DEEPER into the solid and never brackets, silently
                // degrading every such landing to the exact-snap
                // fallback. Hence: (1) the search direction is
                // -sign(phi0)*grad/|grad| -- toward the solid from
                // fluid, toward fluid from solid, always toward the
                // root; (2) the bracket step is seeded at a GEOMETRIC
                // scale guaranteed to cross for a field with
                // |grad phi_r| ~ 1, |phi0| + r (r = this point's own
                // smoothing radius, the worst-case magnitude by
                // which phi_r can differ from the exact field, per
                // Geometry.hpp's derivation) rather than |phi0| alone,
                // with doubling (capped at 12, a safety margin over the
                // single-shot estimate, not the primary bracket
                // mechanism) as a secondary net.
                std::vector<Vec3> descDir(nc, Vec3{0, 0, 0});
                std::vector<char> hasDir(nc, 0);
                for (std::size_t k = 0; k < nc; ++k) {
                    const double gnorm = norm(grad0[k]);
                    if (gnorm > 1e-300) {
                        const double sgn = phi0[k] >= 0.0 ? 1.0 : -1.0;
                        descDir[k] = grad0[k] * (-sgn / gnorm);
                        hasDir[k] = 1;
                    }
                }

                std::vector<double> hiS(nc);
                std::vector<double> aVal(nc), bVal(nc); // bracket endpoints in s
                std::vector<double> fa = phi0, fb(nc, 0.0); // phi at the endpoints
                std::vector<char> bracketed(nc, 0);
                std::vector<long> evalCount(nc, 1); // the phi0 evaluation above counts as eval #1
                // REACH: phi_r differs from the exact field by at most r, so the
                // root lies within (exact remaining distance + r) of the start;
                // twice that allows for an oblique descent direction. A sign
                // change further out is some OTHER surface. MEASURED,
                // win_desc_slot: on the V-shaped field of a duplicated CAD face
                // (no sign change at the wall) the doubling below walked 1.01 m
                // through the pier to its far face and "landed" there. Past the
                // reach a candidate is left unbracketed: the counted exact-snap
                // fallback.
                std::vector<double> reach(nc);
                for (std::size_t k = 0; k < nc; ++k) {
                    hiS[k] = std::fabs(phi0[k]) + r0[k] + 1e-12;
                    aVal[k] = 0.0;
                    const double remaining = norm(exactSnap[candIdx[k]] - pts0[k]);
                    reach[k] = 2.0 * (std::max(std::fabs(phi0[k]), remaining) + r0[k]) + 1e-12;
                }
                constexpr int kMaxBracketDoublings = 12;
                for (int tries = 0; tries < kMaxBracketDoublings; ++tries) {
                    std::vector<std::size_t> active;
                    for (std::size_t k = 0; k < nc; ++k) {
                        if (hasDir[k] && !bracketed[k] && hiS[k] <= reach[k]) active.push_back(k);
                    }
                    if (active.empty()) break;
                    std::vector<Vec3> qp(active.size());
                    std::vector<double> qr(active.size());
                    for (std::size_t a = 0; a < active.size(); ++a) {
                        qp[a] = pts0[active[a]] + descDir[active[a]] * hiS[active[a]];
                        qr[a] = r0[active[a]];
                    }
                    const std::vector<double> qphi = smoothedSignedDistanceBatch(
                        bins, perStlTris[static_cast<std::size_t>(spec.stlIndex)], locationInMesh, qp, qr,
                        smoothCfg);
                    for (std::size_t a = 0; a < active.size(); ++a) {
                        const std::size_t k = active[a];
                        ++evalCount[k];
                        if ((fa[k] > 0.0) != (qphi[a] > 0.0)) {
                            bVal[k] = hiS[k];
                            fb[k] = qphi[a];
                            bracketed[k] = 1;
                        } else {
                            hiS[k] *= 2.0;
                        }
                    }
                }

                // Illinois regula falsi, capped at kMaxRegulaFalsiIters
                // as a BACKSTOP only. Convergence test: |phi(c)| below an
                // absolute tolerance scaled by the wall layer's own
                // thickness (tStep) -- NOT bracket-width shrinkage (the
                // bisection convention CutData.cpp uses): Illinois lets
                // one endpoint's POSITION stick for many iterations even
                // after the ROOT ESTIMATE c has effectively converged
                // (only that endpoint's stored function VALUE gets
                // halved, not its position), so gating on bracket width
                // would keep iterating long after the answer stopped
                // moving -- measured: with a bracket-width gate this
                // averaged ~18 evaluations/candidate, defeating the
                // whole point of the fix.
                constexpr double kLandingAbsTol = 1e-9;
                constexpr int kMaxRegulaFalsiIters = 25;
                std::vector<int> side(nc, 0); // 0 = none yet, 1 = a retained last, 2 = b retained last
                std::vector<double> rootS(nc, 0.0);
                std::vector<char> rfConverged(nc, 0);
                for (std::size_t k = 0; k < nc; ++k) {
                    // Seed the estimate from the initial bracket, in case
                    // it is already converged (loop below never touches
                    // this point).
                    if (hasDir[k] && bracketed[k]) {
                        rootS[k] = bVal[k] - fb[k] * (bVal[k] - aVal[k]) / (fb[k] - fa[k]);
                    }
                }
                const double landTol = kLandingAbsTol * std::max(tStep, 1e-12);
                for (int it = 0; it < kMaxRegulaFalsiIters; ++it) {
                    std::vector<std::size_t> active;
                    for (std::size_t k = 0; k < nc; ++k) {
                        if (hasDir[k] && bracketed[k] && !rfConverged[k]) active.push_back(k);
                    }
                    if (active.empty()) break;
                    std::vector<Vec3> qp(active.size());
                    std::vector<double> qr(active.size());
                    std::vector<double> cVal(active.size());
                    for (std::size_t a = 0; a < active.size(); ++a) {
                        const std::size_t k = active[a];
                        // False-position estimate (linear interpolation
                        // of the root between the two bracket samples).
                        cVal[a] = bVal[k] - fb[k] * (bVal[k] - aVal[k]) / (fb[k] - fa[k]);
                        rootS[k] = cVal[a];
                        qp[a] = pts0[k] + descDir[k] * cVal[a];
                        qr[a] = r0[k];
                    }
                    const std::vector<double> qphi = smoothedSignedDistanceBatch(
                        bins, perStlTris[static_cast<std::size_t>(spec.stlIndex)], locationInMesh, qp, qr,
                        smoothCfg);
                    for (std::size_t a = 0; a < active.size(); ++a) {
                        const std::size_t k = active[a];
                        ++evalCount[k];
                        const double fc = qphi[a];
                        if (std::fabs(fc) < landTol * pointScale[candIdx[k]]) {
                            rfConverged[k] = 1;
                            continue;
                        }
                        if ((fc > 0.0) == (fa[k] > 0.0)) {
                            // `a` side replaced -- if `b` was ALSO the
                            // retained side last time, it is sticking:
                            // halve its stored value (Illinois safeguard)
                            // so the next interpolation weights toward
                            // the other side instead of stalling there.
                            if (side[k] == 2) fb[k] *= 0.5;
                            aVal[k] = cVal[a];
                            fa[k] = fc;
                            side[k] = 1;
                        } else {
                            if (side[k] == 1) fa[k] *= 0.5;
                            bVal[k] = cVal[a];
                            fb[k] = fc;
                            side[k] = 2;
                        }
                    }
                }

                long stepEvalSum = 0, stepEvalMax = 0, stepCandCount = 0;
                for (std::size_t k = 0; k < nc; ++k) {
                    const std::size_t i = candIdx[k];
                    if (hasDir[k] && bracketed[k]) {
                        finalPos[i] = pts0[k] + descDir[k] * rootS[k];
                        stepEvalSum += evalCount[k];
                        stepEvalMax = std::max(stepEvalMax, evalCount[k]);
                        ++stepCandCount;
                    } else {
                        // Fell back to the exact snap -- COUNT it, and
                        // COUNT WHICH of the two distinct causes it was
                        // (a silent fallback, or an unattributed one,
                        // would reintroduce the landing error invisibly):
                        // the gradient estimate was degenerate at the
                        // pre-landing sample itself (never had a search
                        // direction), vs. a direction existed but
                        // kMaxBracketDoublings ran out without a sign
                        // change (the field genuinely did not cross zero
                        // within the geometric-scale step + doublings
                        // tried -- e.g. a very short/thin feature).
                        finalPos[i] = exactSnap[i];
                        if (!hasDir[k]) {
                            ++stats.smoothLandingDegenerateGradFallback;
                        } else {
                            ++stats.smoothLandingNoBracketFallback;
                        }
                        ++stats.smoothLandingBisectFallback; // == the two counters' sum, kept for compatibility
                    }
                }
                stats.smoothLandingEvalSum += stepEvalSum;
                stats.smoothLandingEvalMax = std::max(stats.smoothLandingEvalMax, stepEvalMax);
                stats.smoothLandingCandidateCount += stepCandCount;
                // Every landing ATTEMPT this pass/step is either a
                // counted success (added to smoothLandingCandidateCount
                // above) or a counted fallback -- so
                // smoothLandingCandidateCount + smoothLandingBisectFallback
                // is the comparable denominator for a fallback RATE (both
                // sides accumulate the SAME way, once per attempt, across
                // every quality-gate pass the march re-runs;
                // neither is a per-run total by itself).
            }
        }

        clk.lap("landing");
        // --- "What smoothing destroys must be loud". Only meaningful (and only paid for) on the
        // WALL layer's landed positions, and only when smoothingOn -- this
        // measures the landing above, and (since the seal guard below)
        // also FEEDS it: a face flagged `sealed` is not emitted.
        //
        // Declared outside the block so the validation loop can read it
        // (empty/zero on every step that does not evaluate it, which is
        // exactly the steps where no landing happened).
        std::vector<char> faceSealed(static_cast<std::size_t>(nTop), 0);
        // Sealed faces whose TOP (where a face refused at the landing stays) is
        // itself buried in the solid -- see the seal guard below.
        std::vector<char> faceBuriedTop(static_cast<std::size_t>(nTop), 0);
        if (smoothingOn && lastStep) {
            std::vector<Vec3> landedPts(nFront);
            for (std::size_t i = 0; i < nFront; ++i) landedPts[i] = finalPos[i];

            // Achieved stack thickness: origin's step-0
            // top position (the exact d=t cut position -- `out.points`
            // is append-only, so the ORIGIN point index still holds it
            // unmodified) to this stack's final landed bottom. The cut
            // stays at exact d=t (CutData.cpp untouched) while the wall
            // now moves to phi_r=0, so this measures the confound
            // directly rather than assuming its size.
            for (std::size_t i = 0; i < nFront; ++i) {
                const int origin = originOfPoint(frontPoints[i]);
                const Vec3& topPos = out.points[static_cast<std::size_t>(origin)];
                const double achieved = norm(topPos - landedPts[i]);
                stats.stackThickness.push_back(achieved);
                if (achieved < 0.5 * t * pointScale[i]) ++stats.stackThicknessBelowHalfT;
                if (achieved < 0.25 * t * pointScale[i]) ++stats.stackThicknessBelowQuarterT;
            }

            // Exact signed distance at the landed position: one more
            // closestPointOnSoup per point (the "expensive" exact query
            // the smoothed march otherwise avoids) plus ONE
            // batched classifyVertices call for the sign -- paid only
            // here, once per wall face, not per march step.
            std::vector<double> phiExact(nFront);
            std::vector<Vec3> landedClosest(nFront);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 512)
#endif
            for (std::ptrdiff_t fi = 0; fi < static_cast<std::ptrdiff_t>(nFront); ++fi) {
                const std::size_t i = static_cast<std::size_t>(fi);
                const ClosestHit hit = closestPointOnSoup(bins, landedPts[i]);
                phiExact[i] = std::sqrt(hit.distSq);
                landedClosest[i] = hit.point;
            }
            const std::vector<bool> landedSolid =
                classifyVertices(landedPts, perStlTris[static_cast<std::size_t>(spec.stlIndex)], locationInMesh, &bins);
            for (std::size_t i = 0; i < nFront; ++i) {
                if (landedSolid[i]) phiExact[i] = -phiExact[i];
            }
            const std::vector<double> phiR = smoothedSignedDistanceBatch(
                bins, perStlTris[static_cast<std::size_t>(spec.stlIndex)], locationInMesh, landedPts, radii,
                smoothCfg);

            std::vector<int> ufParent(nFront);
            for (std::size_t i = 0; i < nFront; ++i) ufParent[i] = static_cast<int>(i);
            std::function<int(int)> ufFind = [&](int x) {
                while (ufParent[static_cast<std::size_t>(x)] != x) {
                    ufParent[static_cast<std::size_t>(x)] = ufParent[static_cast<std::size_t>(ufParent[static_cast<std::size_t>(x)])];
                    x = ufParent[static_cast<std::size_t>(x)];
                }
                return x;
            };
            auto ufUnion = [&](int a, int b) {
                a = ufFind(a);
                b = ufFind(b);
                if (a != b) ufParent[static_cast<std::size_t>(a)] = b;
            };
            for (int fi = 0; fi < nTop; ++fi) {
                const IntSpan fp = wallBucket.pointsOf(fi);
                const int lvl = levelOfCell(wallBucket.owner[static_cast<std::size_t>(fi)]);
                const double hLocal = dx0 / static_cast<double>(1 << lvl);
                std::vector<Vec3> loop;
                loop.reserve(static_cast<std::size_t>(fp.size()));
                bool moved = false, sealed = false, buried = false;
                double faceMaxDiff = 0.0;
                int prevIdx = -1, firstIdx = -1;
                for (int p : fp) {
                    loop.push_back(out.points[static_cast<std::size_t>(p)]);
                    const int idx = frontIdx[p];
                    const double diff = std::fabs(phiR[static_cast<std::size_t>(idx)] - phiExact[static_cast<std::size_t>(idx)]);
                    faceMaxDiff = std::max(faceMaxDiff, diff);
                    if (diff > 0.1 * hLocal) moved = true;
                    // The sign disagreement must be BACKED BY DEPTH.
                    // A landed wall face sits ON the surface by construction, so
                    // phiExact there is ~0 and its SIGN is numerical noise -- an
                    // ungated sign test fires on ordinary healthy landings,
                    // swamps the real signal, and grossly overstates sealing,
                    // making the metric useless for the radius decision it
                    // exists to inform.
                    //
                    // Gate on |phiExact| -- the DEPTH the smoothed surface put
                    // this landed face inside the solid -- at the same 0.1*h
                    // `moved` uses, so there is one number to reason about. This
                    // is the stricter choice: with phiExact < 0 < phiR,
                    // |phiR - phiExact| >= |phiExact| always, so gating on the
                    // difference instead would flag strictly more faces.
                    //
                    // A point buried on the side it came from is ACCEPTED: the
                    // smoothed zero rounds a sharp convex edge and lies inside
                    // the solid there, which is the landing working as
                    // designed, only deeper than 0.1*h (MEASURED on a
                    // hydrofoil: 9095 leading/trailing-edge faces, at most
                    // 0.2 h deep; on win_desc_slot 48 faces of a plate
                    // thinner than the smoothing diameter, 0.15 h; refusing
                    // them left each stack a layer short). It is reported as
                    // a buried landing. Refused: a point deeper than
                    // kMaxBuriedDepthFrac * h, a point whose nearest surface
                    // is no longer the one it approached -- it passed the
                    // midline of a rib, where the stack landing from the far
                    // side could cross it -- and, as before, a sign flip on
                    // the fluid side.
                    const double phiE = phiExact[static_cast<std::size_t>(idx)];
                    if ((phiR[static_cast<std::size_t>(idx)] > 0.0) != (phiE > 0.0) &&
                        std::fabs(phiE) > 0.1 * hLocal) {
                        const Vec3& landed = landedPts[static_cast<std::size_t>(idx)];
                        const Vec3& from = out.points[static_cast<std::size_t>(frontPoints[static_cast<std::size_t>(idx)])];
                        if (phiE < 0.0 && -phiE <= kMaxBuriedDepthFrac * hLocal &&
                            dot(landedClosest[static_cast<std::size_t>(idx)] - landed, from - landed) > 0.0) {
                            buried = true;
                            stats.buriedLandingMaxDepth = std::max(stats.buriedLandingMaxDepth, -phiE);
                            stats.buriedLandingMaxDepthOverH = std::max(stats.buriedLandingMaxDepthOverH, -phiE / hLocal);
                        } else {
                            sealed = true;
                        }
                    }
                    if (firstIdx < 0) firstIdx = idx;
                    if (prevIdx >= 0) ufUnion(prevIdx, idx);
                    prevIdx = idx;
                }
                if (firstIdx >= 0 && prevIdx >= 0) ufUnion(prevIdx, firstIdx); // close the loop
                const double area = loopArea(loop);
                stats.smoothResidAreaSum += area;
                // Reported in units of h_local: each face's
                // own hLocal normalizes its own contribution BEFORE the
                // area weighting, so smoothResidWeightedSum /
                // smoothResidAreaSum is directly the area-weighted mean
                // of |phi_r - phi| / h_local across the whole wall.
                stats.smoothResidWeightedSum += area * (faceMaxDiff / hLocal);
                stats.smoothResidMaxOverH = std::max(stats.smoothResidMaxOverH, faceMaxDiff / hLocal);
                if (moved) {
                    ++stats.smoothMovedFaces;
                    stats.smoothMovedArea += area;
                }
                faceSealed[static_cast<std::size_t>(fi)] = sealed ? 1 : 0;
                if (buried && !sealed) {
                    ++stats.buriedLandingFaces;
                    stats.buriedLandingArea += area;
                }
            }
            // A sealed face is refused at the LANDING only: it reverts to a wall
            // face at the interface it had already reached (its top loop). That
            // is sound only if the top is not itself inside the solid, so the
            // same depth predicate is applied there -- to the sealed faces'
            // points only, one exact query + one batched parity pass. A buried
            // top falls back to condemning the whole stack, and only those
            // count as sealed regions: the metric reports what SHIPS.
            {
                std::vector<int> sealedIdx;
                std::unordered_map<int, int> slotOf;
                for (int fi = 0; fi < nTop; ++fi) {
                    if (!faceSealed[static_cast<std::size_t>(fi)]) continue;
                    for (int p : wallBucket.pointsOf(fi)) {
                        const int ix = frontIdx[p];
                        if (slotOf.emplace(ix, static_cast<int>(sealedIdx.size())).second) sealedIdx.push_back(ix);
                    }
                }
                std::vector<Vec3> topPts(sealedIdx.size());
                std::vector<double> phiTop(sealedIdx.size());
                for (std::size_t k = 0; k < sealedIdx.size(); ++k) {
                    topPts[k] = out.points[static_cast<std::size_t>(frontPoints[static_cast<std::size_t>(sealedIdx[k])])];
                    phiTop[k] = std::sqrt(closestPointOnSoup(bins, topPts[k]).distSq);
                }
                const std::vector<bool> topSolid =
                    classifyVertices(topPts, perStlTris[static_cast<std::size_t>(spec.stlIndex)], locationInMesh, &bins);
                for (int fi = 0; fi < nTop; ++fi) {
                    if (!faceSealed[static_cast<std::size_t>(fi)]) continue;
                    const double hLocal = dx0 / static_cast<double>(1 << levelOfCell(wallBucket.owner[static_cast<std::size_t>(fi)]));
                    for (int p : wallBucket.pointsOf(fi)) {
                        const std::size_t k = static_cast<std::size_t>(slotOf.at(frontIdx[p]));
                        if (topSolid[k] && phiTop[k] > 0.1 * hLocal) faceBuriedTop[static_cast<std::size_t>(fi)] = 1;
                    }
                }
            }
            // Connected-region accounting over sealed faces only (union-
            // find already unioned every sealed face's own point loop
            // above -- two sealed faces sharing an edge share two
            // points, hence the same root, by construction).
            std::unordered_map<int, std::size_t> rootToRegion;
            for (int fi = 0; fi < nTop; ++fi) {
                if (!faceBuriedTop[static_cast<std::size_t>(fi)]) continue;
                const IntSpan fp = wallBucket.pointsOf(fi);
                std::vector<Vec3> loop;
                loop.reserve(static_cast<std::size_t>(fp.size()));
                int anyIdx = -1;
                for (int p : fp) {
                    loop.push_back(out.points[static_cast<std::size_t>(p)]);
                    if (anyIdx < 0) anyIdx = frontIdx[p];
                }
                const int root = ufFind(anyIdx);
                auto it = rootToRegion.find(root);
                const double area = loopArea(loop);
                const Vec3 centroid = loopCentroid(loop);
                if (it == rootToRegion.end()) {
                    LayerStats::SealedRegion region;
                    region.count = 1;
                    region.area = area;
                    region.centroid = centroid * area;
                    rootToRegion.emplace(root, stats.smoothSealedRegions.size());
                    stats.smoothSealedRegions.push_back(region);
                } else {
                    LayerStats::SealedRegion& region = stats.smoothSealedRegions[it->second];
                    region.centroid = region.centroid + centroid * area;
                    region.area += area;
                    ++region.count;
                }
            }
            for (LayerStats::SealedRegion& region : stats.smoothSealedRegions) {
                if (region.area > 1e-300) region.centroid = region.centroid * (1.0 / region.area);
            }

            // Opt-in CSV dump (same shape as NINJA_LAYER_GATE_DUMP),
            // one row per landed wall face flagged moved and/or sealed --
            // patch, origin, world centroid, level, phi, phi_r -- so any
            // smoothing surprise is attributable without a rerun.
            const char* smoothDumpPath = std::getenv("NINJA_SMOOTH_DUMP");
            if (smoothDumpPath) {
                static std::FILE* dumpFile = nullptr;
                static bool dumpHeaderWritten = false;
                if (!dumpFile) {
                    dumpFile = std::fopen(smoothDumpPath, "w");
                    if (dumpFile) {
                        std::fprintf(dumpFile, "patchOrdinal,origin,x,y,z,level,phiExact,phiR,moved,sealed\n");
                        dumpHeaderWritten = true;
                    }
                }
                if (dumpFile && dumpHeaderWritten) {
                    for (int fi = 0; fi < nTop; ++fi) {
                        const IntSpan fp = wallBucket.pointsOf(fi);
                        const int lvl = levelOfCell(wallBucket.owner[static_cast<std::size_t>(fi)]);
                        std::vector<Vec3> loop;
                        loop.reserve(static_cast<std::size_t>(fp.size()));
                        double maxDiff = 0.0;
                        double repPhiExact = 0.0, repPhiR = 0.0;
                        const double hLocal = dx0 / static_cast<double>(1 << lvl);
                        bool moved = false;
                        for (int p : fp) {
                            loop.push_back(out.points[static_cast<std::size_t>(p)]);
                            const int idx = frontIdx[p];
                            const double diff =
                                std::fabs(phiR[static_cast<std::size_t>(idx)] - phiExact[static_cast<std::size_t>(idx)]);
                            if (diff >= maxDiff) {
                                maxDiff = diff;
                                repPhiExact = phiExact[static_cast<std::size_t>(idx)];
                                repPhiR = phiR[static_cast<std::size_t>(idx)];
                            }
                            if (diff > 0.1 * hLocal) moved = true;
                        }
                        const bool sealed = faceSealed[static_cast<std::size_t>(fi)] != 0;
                        if (!moved && !sealed) continue;
                        const Vec3 c = loopCentroid(loop);
                        std::fprintf(dumpFile, "%d,%d,%.9f,%.9f,%.9f,%d,%.9g,%.9g,%d,%d\n", pOrd, faceOrigin[static_cast<std::size_t>(fi)],
                                     c.x, c.y, c.z, lvl, repPhiExact, repPhiR, moved ? 1 : 0, sealed ? 1 : 0);
                    }
                    std::fflush(dumpFile);
                }
            }
        }

        clk.lap("smoothing disclosure");
        // --- New bottom points (deterministic creation order: ascending
        // front-point original index, independent of face iteration order).
        std::vector<int> bottomPointIdx(nFront, -1);
        std::vector<int> bottomPointLevelAcc(nFront, -1); // max over incident faces' owner cell level
        for (std::size_t i = 0; i < nFront; ++i) {
            bottomPointIdx[i] = static_cast<int>(out.points.size());
            out.points.push_back(finalPos[i]);
            pointLevel.push_back(0); // filled below once face owners are known
            // Origin inheritance (bookkeeping only -- no geometry here).
            pointOrigin[bottomPointIdx[i]] = originOfPoint(frontPoints[i]);
        }

        clk.lap("new bottom points");
        // --- Prism validation using FINAL positions + forced
        // drops (layersDebug{forceDropSphere}).
        std::vector<bool> faceValid(static_cast<std::size_t>(nTop), true);
        // Candidate prism centres (every face that passes prismValid), so a
        // side quad can be scored against the prism it will share it with.
        std::vector<char> candOk(static_cast<std::size_t>(nTop), 0);
        std::vector<Vec3> candCentroid(static_cast<std::size_t>(nTop));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256)
#endif
        for (int fi = 0; fi < nTop; ++fi) {
            if (faceHeld(fi)) continue;
            const IntSpan fp = wallBucket.pointsOf(fi);
            std::vector<Vec3> topLoop, botLoop;
            for (int p : fp) {
                topLoop.push_back(out.points[static_cast<std::size_t>(p)]);
                botLoop.push_back(out.points[static_cast<std::size_t>(bottomPointIdx[static_cast<std::size_t>(frontIdx.at(p))])]);
            }
            const double fs = faceScale[static_cast<std::size_t>(fi)];
            if (!prismValid(topLoop, botLoop, minHeight * fs, minVolEps * fs * fs * fs, nullptr, kMaxPrismAspect,
                            minAchievedHeightFrac() * tStep * fs, nullptr)) {
                continue;
            }
            const std::size_t n = topLoop.size();
            std::vector<std::vector<Vec3>> loops;
            loops.push_back(std::vector<Vec3>(topLoop.rbegin(), topLoop.rend()));
            loops.push_back(botLoop);
            for (std::size_t i = 0; i < n; ++i) {
                loops.push_back({topLoop[i], topLoop[(i + 1) % n], botLoop[(i + 1) % n], botLoop[i]});
            }
            candCentroid[static_cast<std::size_t>(fi)] = polyhedronCentroidFromOutwardLoops(loops);
            candOk[static_cast<std::size_t>(fi)] = 1;
        }
        auto sideNeighboursOf = [&](int fi) {
            const IntSpan fp = wallBucket.pointsOf(fi);
            const int n = fp.size();
            std::vector<const Vec3*> nb(static_cast<std::size_t>(n), nullptr);
            for (int i = 0; i < n; ++i) {
                auto it = edgeToFaces.find(makeEdgeKey(fp[i], fp[(i + 1) % n]));
                if (it == edgeToFaces.end() || it->second.size() != 2) continue;
                const int other = it->second[0] == fi ? it->second[1] : it->second[0];
                if (candOk[static_cast<std::size_t>(other)]) nb[static_cast<std::size_t>(i)] = &candCentroid[static_cast<std::size_t>(other)];
            }
            return nb;
        };
        // In-march quality guard (the POST-CONDITION): every prism is
        // evaluated with the gate's own checkMesh-mimicking predicate
        // (face-pyramid orientation + tet decomposition, owner side of
        // the top face included) BEFORE it is committed; a failing prism
        // is never emitted -- its face is reverted (terraced) at this
        // step exactly like a prismValid failure, and the mesh simply
        // stays away from the wall there. MEASURED necessary
        // (bm_layers_wfp): the re-march gate below is capped at
        // kMaxGatePasses and never converged on that terrain (443, 218,
        // 189, 144, 127, 148 newly condemned stacks per pass), so 148
        // known-bad stacks SHIPPED. The gate loop is kept as a safety
        // net; with this guard it should find nothing.
        std::vector<char> qualityBad(static_cast<std::size_t>(nTop), 0);
        std::vector<int> dropReason(static_cast<std::size_t>(nTop), 0); // LayerDropReason; 0 = not refused itself
        std::vector<LayerCellQuality> faceQuality(static_cast<std::size_t>(nTop));
        for (int fi = 0; fi < nTop; ++fi) {
            const IntSpan fp = wallBucket.pointsOf(fi);
            std::vector<Vec3> topLoop, botLoop;
            topLoop.reserve(static_cast<std::size_t>(fp.size()));
            botLoop.reserve(static_cast<std::size_t>(fp.size()));
            for (int p : fp) {
                topLoop.push_back(out.points[static_cast<std::size_t>(p)]);
                botLoop.push_back(out.points[static_cast<std::size_t>(bottomPointIdx[static_cast<std::size_t>(frontIdx[p])])]);
            }
            int fcode = 0;
            double meanHeight = -1.0;
            const double fs = faceScale[static_cast<std::size_t>(fi)];
            const double tStepF = tStep * fs; // this face's own layer step
            bool valid = prismValid(topLoop, botLoop, minHeight * fs, minVolEps * fs * fs * fs, &fcode, kMaxPrismAspect,
                                    minAchievedHeightFrac() * tStepF, &meanHeight);
            // A held face (see faceHeld) is never extruded -- it stays a wall
            // face for the rest of the march. MEASURED (bm_layers_wfp):
            // re-extruding the terrace seams always failed, every step --
            // 19474 -> 8979 dropped faces from holding them.
            // fcode 99 keeps a held face out of every guard's own accounting.
            if (faceHeld(fi)) {
                valid = false;
                fcode = 99;
                meanHeight = -1.0;
            }
            if (!valid) dropReason[static_cast<std::size_t>(fi)] = fcode == 99 ? kDropHeldSeam : fcode;
            if (meanHeight >= 0.0 && tStepF > 0.0) {
                const double frac = meanHeight / tStepF;
                if (layerHeightHist() && valid) {
                    // bin index 0..19 over [0,1), 20 for >= 1
                    const int b = frac >= 1.0 ? 20 : static_cast<int>(frac * 20.0);
                    ++g_heightHistBins[static_cast<std::size_t>(std::min(9, layerJ))]
                                      [static_cast<std::size_t>(std::max(0, std::min(20, b)))];
                }
                // The shipped-mesh number: the worst achieved/requested
                // layer height over every prism the march actually
                // EMITS, on the wall step (layerJ == 0, the first layer
                // -- the one whose collapse pins the solver's Courant
                // number). Gated per case as min_layer_height_frac_above.
                if (valid && lastStep) {
                    stats.minLayerHeightFrac = std::min(stats.minLayerHeightFrac, frac);
                }
            }
            if (!valid && fcode == 10) {
                ++stats.thinDroppedFaces;
                stats.thinDroppedArea += loopArea(topLoop);
                // TERRACE, do not condemn the stack. MEASURED:
                // condemning the whole stack the way
                // the seal guard does cascades -- removing a stack
                // perturbs its neighbours' march, which produces more
                // short prisms, which condemn more stacks: 372 / 651 /
                // 201 / 75 / 29 / 10 over six passes, the gate's cap
                // reached, 1216 stacks and 19.3 m2 of wall gone on a
                // window of 12 m2 of healthy layers. Reverting only the
                // FAILING STEP -- the ordinary prismValid treatment,
                // which failCodes 2, 8 and 9 already take -- removes the
                // collapsed cell and nothing else: that face simply
                // gets fewer, thicker layers and the wall stays where
                // the march had already reached. Graceful degradation
                // is the design's own answer here and it is the cheaper
                // one; whole-stack condemnation is reserved for the
                // seal guard, where the defect is the landing itself
                // and no partial retreat is meaningful.
                const int originT = faceGateStack[static_cast<std::size_t>(fi)];
                if (originT >= 0) {
                    report.thinStacks.insert({pOrd, originT});
                } else {
                    ++report.orphanThinFaces;
                }
            }
            // --- Seal guard (the FACING-WALL post-condition). A landed
            // wall face whose landed points sit INSIDE the solid by a
            // geometrically meaningful depth is the signature of the
            // smoothed offset surface having crossed the wall: in a gap
            // narrower than about twice the smoothing radius, phi_r's
            // zero level set is displaced into the solid on BOTH facing
            // walls, the two landings cross, and the passage is sealed
            // by prisms. The test is generic (it is the sign/depth
            // disagreement between the smoothed and the exact field at
            // the landed position, not a named feature), it is the SAME
            // predicate the sealed-region instrument reports, and it
            // costs nothing new -- both fields are already evaluated
            // there.
            //
            // The action is the contract's: never seal. The stack is
            // dropped (the face reverts to the offset cut and is
            // condemned below, so the gate's next pass never extrudes it
            // at all) -- the mesh simply stays away from the wall and
            // the passage keeps whatever the exact d=t offset cut left
            // of it. Shrinking instead of dropping was rejected: the
            // depth by which phi_r crossed is not a step length the
            // march can back off along (the whole stack marched on the
            // crossed field), and a partial retreat leaves a prism whose
            // bottom is still nearer the OPPOSITE wall than to its own.
            if (valid && faceSealed[static_cast<std::size_t>(fi)]) {
                valid = false;
                dropReason[static_cast<std::size_t>(fi)] = kDropSealed;
                ++stats.sealDroppedFaces;
                const IntSpan fpS = wallBucket.pointsOf(fi);
                std::vector<Vec3> loopS;
                loopS.reserve(static_cast<std::size_t>(fpS.size()));
                for (int p : fpS) loopS.push_back(out.points[static_cast<std::size_t>(p)]);
                stats.sealDroppedArea += loopArea(loopS);
                const int originS = faceGateStack[static_cast<std::size_t>(fi)];
                if (originS >= 0) {
                    // Refuse only the LANDING step -- the face reverts to a
                    // wall face at the last interface it reached, exactly as
                    // the collapsed-prism guard does. The whole stack goes
                    // back to the offset cut only when that interface is
                    // itself buried. MEASURED (bm_layers_wfp): condemning
                    // every sealed stack re-marched with a full-height hole
                    // per stack and doubled the dropped faces around each.
                    if (faceBuriedTop[static_cast<std::size_t>(fi)]) report.badStacks.insert({pOrd, originS});
                    report.sealedStacks.insert({pOrd, originS});
                } else {
                    ++report.orphanSealedFaces;
                }
            }
            if (valid) {
                const Vec3 c = loopCentroid(topLoop);
                for (const ForceDropSphere& sph : forceDropSpheres) {
                    if (norm(c - sph.centre) <= sph.radius) {
                        valid = false;
                        dropReason[static_cast<std::size_t>(fi)] = kDropForced;
                        break;
                    }
                }
            }
            if (valid) {
                const int coreCell = wallBucket.owner[static_cast<std::size_t>(fi)];
                const Vec3* ownerCentroid = lqOwnerCentroidOf(coreCell);
                const std::vector<const Vec3*> sideNb = sideNeighboursOf(fi);
                const LayerCellQuality q = evaluateLayerCellQuality(topLoop, botLoop,
                                                                    tStep * faceScale[static_cast<std::size_t>(fi)],
                                                                    ownerCentroid, &sideNb);
                faceQuality[static_cast<std::size_t>(fi)] = q;
                if (q.nBadPyramidFaces > 0 || q.nBadTetFaces > 0 || q.nHighSkewFaces > 0 || q.nHighNonOrthFaces > 0) {
                    qualityBad[static_cast<std::size_t>(fi)] = 1;
                    dropReason[static_cast<std::size_t>(fi)] = q.nBadPyramidFaces > 0 ? kDropPyramid
                                                               : q.nBadTetFaces > 0   ? kDropTet
                                                               : q.nHighSkewFaces > 0 ? kDropSkew
                                                                                      : kDropNonOrth;
                    ++stats.qualityDroppedFaces;
                    if (gateDump) {
                        const int origin = faceGateStack[static_cast<std::size_t>(fi)];
                        report.condemnedSites.push_back(
                            {pOrd, origin, q.centroid.x, q.centroid.y, q.centroid.z, levelOfCell(coreCell)});
                    }
                }
            }
            faceValid[static_cast<std::size_t>(fi)] = valid;
        }
        // The reverted set: every failing face, and only those. Its healthy
        // neighbours are NOT dropped with it -- the seam quads close a
        // reverted face against a full-height neighbour exactly as they do
        // for a gate removal. MEASURED: a one-ring dilation here dropped
        // 46% of all faces on a hydrofoil window (20100 -> 10769 without
        // it) and halved the drops on every layered benchmark, while the
        // fluid volume moved TOWARD the layers-off reference
        // (win_wfp_coarse_layers 27.07 -> 30.54, reference 30.77).
        std::vector<bool> reverted(static_cast<std::size_t>(nTop));
        for (int fi = 0; fi < nTop; ++fi) reverted[static_cast<std::size_t>(fi)] = !faceValid[static_cast<std::size_t>(fi)];
        // Gate removals: stacks condemned by an earlier pass are
        // reverted on EVERY step, so they never extrude at all -- same
        // code path (and therefore the same seam closure, the same
        // conformal terracing against full-height neighbours, and the
        // same "infeasible wall area" accounting) as a face the march
        // drops on its own: a condemned stack's neighbours are removed
        // only if they in turn produce a bad cell, which the next pass
        // measures. What the seams DO need is a gate identity, see
        // faceGateStack.
        if (!gateRemoved.empty()) {
            for (int fi = 0; fi < nTop; ++fi) {
                const int o = faceGateStack[static_cast<std::size_t>(fi)];
                if (o >= 0 && gateRemoved.count({pOrd, o})) reverted[static_cast<std::size_t>(fi)] = true;
            }
        }
        for (int fi = 0; fi < nTop; ++fi) {
            if (qualityBad[static_cast<std::size_t>(fi)]) reverted[static_cast<std::size_t>(fi)] = true;
        }
        // A side quad scored against its neighbour's prism becomes a SEAM
        // (boundary) face when that neighbour is reverted after all, and a
        // boundary face is held to its own prism's centre alone. Re-score
        // those until no new face reverts (MEASURED, bm_layers_wfp: one
        // seam quad at skewness 4.7 shipped without this).
        for (bool changed = true; changed;) {
            changed = false;
            for (int fi = 0; fi < nTop; ++fi) {
                if (reverted[static_cast<std::size_t>(fi)]) continue;
                const IntSpan fp = wallBucket.pointsOf(fi);
                const int n = fp.size();
                const Vec3& cc = faceQuality[static_cast<std::size_t>(fi)].centroid;
                for (int i = 0; i < n; ++i) {
                    auto it = edgeToFaces.find(makeEdgeKey(fp[i], fp[(i + 1) % n]));
                    if (it == edgeToFaces.end() || it->second.size() != 2) continue;
                    const int other = it->second[0] == fi ? it->second[1] : it->second[0];
                    if (!candOk[static_cast<std::size_t>(other)] || !reverted[static_cast<std::size_t>(other)]) continue;
                    const int a = fp[i];
                    const int b = fp[(i + 1) % n];
                    const std::vector<Vec3> quad{
                        out.points[static_cast<std::size_t>(a)], out.points[static_cast<std::size_t>(b)],
                        out.points[static_cast<std::size_t>(bottomPointIdx[static_cast<std::size_t>(frontIdx.at(b))])],
                        out.points[static_cast<std::size_t>(bottomPointIdx[static_cast<std::size_t>(frontIdx.at(a))])]};
                    const bool skewBad = faceSkewnessOf(quad, cc, nullptr) > kSkewMax;
                    if (skewBad || faceNonOrthDegOf(quad, cc, nullptr) > kNonOrthMaxDeg) {
                        reverted[static_cast<std::size_t>(fi)] = true;
                        qualityBad[static_cast<std::size_t>(fi)] = 1;
                        dropReason[static_cast<std::size_t>(fi)] = skewBad ? kDropSkew : kDropNonOrth;
                        ++stats.qualityDroppedFaces;
                        changed = true;
                        break;
                    }
                }
            }
        }
        // --- A prism that fails below the first layer is merged into the
        // prism above it instead of being reverted: that cell grows down to
        // the new bottom and the interface between the two, the face whose
        // warp the thin prism could not carry, is not emitted. The stack
        // keeps one layer fewer but still reaches the wall. MEASURED: the
        // wall-layer prisms refused on their top face's tets sit under an
        // interface warped by 0.6-1.1x their own height; no triangulation
        // helps (the prism centre stays behind part of the face), while the
        // merged cell passed every check for 172 of 177 on bm_layers_wfp.
        // The merged cell is held to the same tests as any prism: tets and
        // skewness on all its faces, the interface above two-sided, and the
        // side faces shared with a neighbour two-sided (skewness and
        // non-orthogonality) against the FINAL neighbours -- iterated, since
        // a neighbour may merge too.
        std::vector<int> mergeInto(static_cast<std::size_t>(nTop), -1);
        std::vector<Vec3> mergedCentroid(static_cast<std::size_t>(nTop));
        {
            const auto P = [&](int p) -> const Vec3& { return out.points[static_cast<std::size_t>(p)]; };
            const auto botOf = [&](int p) { return bottomPointIdx[static_cast<std::size_t>(frontIdx.at(p))]; };
            // Upper side i spans T -> M, lower side i spans M -> bottom.
            struct MergedCell {
                std::vector<Vec3> top;                     // parent's top loop, stored orientation
                std::vector<std::vector<Vec3>> loops;      // outward: top reversed, bottom, then upper/lower per side
            };
            const auto buildMerged = [&](int fi, const std::vector<int>& T) {
                MergedCell mc;
                const IntSpan M = wallBucket.pointsOf(fi);
                const std::size_t n = T.size();
                for (int p : T) mc.top.push_back(P(p));
                mc.loops.push_back(std::vector<Vec3>(mc.top.rbegin(), mc.top.rend()));
                std::vector<Vec3> bot;
                for (int p : M) bot.push_back(P(botOf(p)));
                mc.loops.push_back(bot);
                for (std::size_t i = 0; i < n; ++i) {
                    const std::size_t j = (i + 1) % n;
                    mc.loops.push_back({P(T[i]), P(T[j]), P(M[j]), P(M[i])});
                    mc.loops.push_back({P(M[i]), P(M[j]), P(botOf(M[j])), P(botOf(M[i]))});
                }
                return mc;
            };
            const auto eligible = [&](int fi) -> const std::vector<int>* {
                const int r = dropReason[static_cast<std::size_t>(fi)];
                if (r == 0 || r == kDropHeldSeam || r == kDropSealed || r == kDropForced) return nullptr;
                const int gs = faceGateStack[static_cast<std::size_t>(fi)];
                if (gs >= 0 && gateRemoved.count({pOrd, gs})) return nullptr;
                const int owner = wallBucket.owner[static_cast<std::size_t>(fi)];
                if (mergedPrisms.count(owner)) return nullptr;
                const auto itT = prismTopOf.find(owner);
                const auto itB = prismBottomOf.find(owner);
                if (itT == prismTopOf.end() || itB == prismBottomOf.end()) return nullptr;
                const IntSpan M = wallBucket.pointsOf(fi);
                // fi must be that prism's own bottom face, not a seam it owns.
                if (itB->second.size() != static_cast<std::size_t>(M.size()) ||
                    !std::equal(itB->second.begin(), itB->second.end(), M.b)) {
                    return nullptr;
                }
                return &itT->second;
            };
            // The checks that do not depend on the neighbours' decisions.
            const auto ownChecks = [&](int fi, const MergedCell& mc, const Vec3& c) {
                LayerCellQuality r;
                const int owner = wallBucket.owner[static_cast<std::size_t>(fi)];
                const auto par = prismParentOf.find(owner);
                const Vec3* parentC = par != prismParentOf.end() ? lqOwnerCentroidOf(par->second) : nullptr;
                if (!parentC) return false;
                if (evalFaceAgainstApex(mc.loops[0], c, r) || evalFaceAgainstApex(mc.top, *parentC, r) ||
                    !faceHasUsableBasePoint(*parentC, &c, mc.top, kMinTetQuality) ||
                    faceSkewnessOf(mc.top, *parentC, &c) > kSkewMax || faceNonOrthDegOf(mc.top, *parentC, &c) > kNonOrthMaxDeg) {
                    return false;
                }
                for (std::size_t k = 1; k < mc.loops.size(); ++k) {
                    if (evalFaceAgainstApex(mc.loops[k], c, r) || !faceHasUsableBasePoint(c, nullptr, mc.loops[k], kMinTetQuality) ||
                        faceSkewnessOf(mc.loops[k], c, nullptr) > kSkewMax) {
                        return false;
                    }
                }
                if (r.nBadPyramidFaces > 0) return false;
                // The part added below must not be folded on its own.
                const double fs = faceScale[static_cast<std::size_t>(fi)];
                const std::size_t n = mc.top.size();
                std::vector<Vec3> mid, bot;
                for (std::size_t i = 0; i < n; ++i) {
                    mid.push_back(mc.loops[2 + 2 * i][3]);
                    bot.push_back(mc.loops[1][i]);
                }
                return prismValid(mid, bot, 0.0, minVolEps * fs * fs * fs);
            };
            // The side faces shared with a neighbour, against the neighbour's
            // final centre (merged, newly built, or unchanged). One-sided
            // skewness does not bound the two-sided one. MEASURED
            // (bm_layers_wfp at grid x8): a merged cell's side passed
            // one-sided and shipped at skewness 4.07 against its neighbour.
            const auto sideBad = [&](const std::vector<Vec3>& loop, const Vec3& c, const Vec3* nc) {
                return nc && (faceNonOrthDegOf(loop, c, nc) > kNonOrthMaxDeg || faceSkewnessOf(loop, c, nc) > kSkewMax);
            };
            const auto sidesOk = [&](int fi, const MergedCell& mc, const Vec3& c) {
                const IntSpan M = wallBucket.pointsOf(fi);
                const int n = M.size();
                const int owner = wallBucket.owner[static_cast<std::size_t>(fi)];
                for (int i = 0; i < n; ++i) {
                    const auto it = edgeToFaces.find(makeEdgeKey(M[i], M[(i + 1) % n]));
                    if (it == edgeToFaces.end()) continue;
                    for (int g : it->second) {
                        if (g == fi) continue;
                        const int nbrCell = wallBucket.owner[static_cast<std::size_t>(g)];
                        const bool gMerged = mergeInto[static_cast<std::size_t>(g)] >= 0;
                        // Upper side: shared with the neighbouring prism of the step above.
                        if (nbrCell != owner) {
                            const Vec3* nc = gMerged ? &mergedCentroid[static_cast<std::size_t>(g)] : lqOwnerCentroidOf(nbrCell);
                            if (sideBad(mc.loops[static_cast<std::size_t>(2 + 2 * i)], c, nc)) return false;
                        }
                        // Lower side: shared with whatever g becomes this step.
                        const Vec3* lc = gMerged ? &mergedCentroid[static_cast<std::size_t>(g)]
                                         : !reverted[static_cast<std::size_t>(g)] ? &faceQuality[static_cast<std::size_t>(g)].centroid
                                                                                  : nullptr;
                        if (sideBad(mc.loops[static_cast<std::size_t>(3 + 2 * i)], c, lc)) return false;
                    }
                }
                return true;
            };
            std::vector<MergedCell> cells(static_cast<std::size_t>(nTop));
            std::vector<int> tentative;
            for (int fi = 0; fi < nTop; ++fi) {
                if (!reverted[static_cast<std::size_t>(fi)]) continue;
                const std::vector<int>* T = eligible(fi);
                if (!T) continue;
                MergedCell mc = buildMerged(fi, *T);
                const Vec3 c = polyhedronCentroidFromOutwardLoops(mc.loops);
                if (!ownChecks(fi, mc, c)) continue;
                mergeInto[static_cast<std::size_t>(fi)] = wallBucket.owner[static_cast<std::size_t>(fi)];
                mergedCentroid[static_cast<std::size_t>(fi)] = c;
                cells[static_cast<std::size_t>(fi)] = std::move(mc);
                tentative.push_back(fi);
            }
            for (bool changed = true; changed;) {
                changed = false;
                for (int fi : tentative) {
                    if (mergeInto[static_cast<std::size_t>(fi)] < 0) continue;
                    if (!sidesOk(fi, cells[static_cast<std::size_t>(fi)], mergedCentroid[static_cast<std::size_t>(fi)])) {
                        mergeInto[static_cast<std::size_t>(fi)] = -1;
                        changed = true;
                    }
                }
            }
            for (int fi : tentative) {
                if (mergeInto[static_cast<std::size_t>(fi)] < 0) continue;
                // No longer a drop: take back what the validation loop recorded.
                if (dropReason[static_cast<std::size_t>(fi)] == kDropCollapsed) {
                    std::vector<Vec3> topL;
                    for (int p : wallBucket.pointsOf(fi)) topL.push_back(P(p));
                    --stats.thinDroppedFaces;
                    stats.thinDroppedArea -= loopArea(topL);
                    report.thinStacks.erase({pOrd, faceGateStack[static_cast<std::size_t>(fi)]});
                }
                if (qualityBad[static_cast<std::size_t>(fi)]) --stats.qualityDroppedFaces;
                qualityBad[static_cast<std::size_t>(fi)] = 0;
                dropReason[static_cast<std::size_t>(fi)] = 0;
                reverted[static_cast<std::size_t>(fi)] = false;
                LayerCellQuality q;
                q.centroid = mergedCentroid[static_cast<std::size_t>(fi)];
                faceQuality[static_cast<std::size_t>(fi)] = q;
                const int owner = mergeInto[static_cast<std::size_t>(fi)];
                if (static_cast<std::size_t>(owner) >= cellReached.size() || cellReached[static_cast<std::size_t>(owner)]) {
                    ++stats.mergedLayerFaces; // counted, like a drop, only where the stack ships
                }
            }
        }
        for (int fi = 0; fi < nTop; ++fi) {
            if (!reverted[static_cast<std::size_t>(fi)]) continue;
            const int ownerHere = wallBucket.owner[static_cast<std::size_t>(fi)];
            if (ownerHere >= 0 && static_cast<std::size_t>(ownerHere) < cellReached.size() &&
                !cellReached[static_cast<std::size_t>(ownerHere)]) {
                continue; // never ships
            }
            ++stats.perStepDropped[static_cast<std::size_t>(step)];
            const int o = faceOrigin[static_cast<std::size_t>(fi)];
            if (o >= 0 && droppedOrigins.insert(o).second) {
                ++stats.droppedFaces;
                int reason = dropReason[static_cast<std::size_t>(fi)];
                const int gs = faceGateStack[static_cast<std::size_t>(fi)];
                if (reason == 0 && gs >= 0 && gateRemoved.count({pOrd, gs})) reason = kDropCondemned;
                if (stats.droppedByReason.empty()) stats.droppedByReason.assign(kNumLayerDropReasons, 0);
                ++stats.droppedByReason[static_cast<std::size_t>(std::min(reason, kNumLayerDropReasons - 1))];
                // Layer-feasibility accounting: a face that never completes ANY layer
                // step, and stays reverted to a plain wall face at the
                // offset-cut position, is exactly the intended
                // "infeasible region keeps the offset cut and seals"
                // case -- reusing the march's OWN existing drop
                // accounting (no new mechanism), just also weighting it
                // by area so the trade is a wall-AREA number, not only a
                // face count.
                const IntSpan fp0 = wallBucket.pointsOf(fi);
                std::vector<Vec3> loop0;
                loop0.reserve(static_cast<std::size_t>(fp0.size()));
                for (int p : fp0) loop0.push_back(out.points[static_cast<std::size_t>(p)]);
                stats.infeasibleWallArea += loopArea(loop0);
                if (gateDump) {
                    const Vec3 fc = areaWeightedCentroid(loop0);
                    const int coreCellHere = wallBucket.owner[static_cast<std::size_t>(fi)];
                    report.terracedFaces.push_back({pOrd, o, fc.x, fc.y, fc.z, levelOfCell(coreCellHere), step,
                                                     nSteps, reason});
                }
            }
        }


        clk.lap("validation + drops");
        // --- Emission: build the new wall bucket (dropped faces copied
        // unchanged; kept faces become internal + get a bottom face),
        // side faces (internal between two kept prisms, else boundary),
        // and prism cells.
        FaceStore newWallBucket;
        std::vector<int> newFaceOrigin; // parallel to newWallBucket (bookkeeping only)
        std::vector<int> newFaceGateStack; // parallel to newWallBucket (see faceGateStack)
        std::vector<int> topFaceToPrism(static_cast<std::size_t>(nTop), -1);
        // Pending shared side faces, indexed by their shared top edge
        // at PUSH time. Without this index the second pass below would
        // LINEARLY SCAN all of internalOut once per shared prism edge
        // -- O(nTopEdges x nInternalFaces), uncompletable on large
        // layered meshes. Every face that scan could match is a pending
        // side face pushed here with neighbour -1, and its point set
        // {a, b, bottom(a), bottom(b)} determines the edge (bottom
        // indices are fresh and distinct from all top indices), so
        // candidates for edge (a,b) are exactly this list, in the same
        // creation order the global scan visited them. The second pass
        // applies the IDENTICAL predicate to this list only.
        std::unordered_map<EdgeKey, std::vector<int>, EdgeKeyHash> pendingSideByEdge;
        for (int fi = 0; fi < nTop; ++fi) {
            // Capture the ORIGINAL wall seed-face
            // non-planarity once, at step 0, for every wall face (kept
            // or dropped) -- step 0's wallBucket IS the un-marched
            // cutter geometry (see faceOrigin's minting above), and a
            // dropped-then-later-flagged origin still needs this value.
            if (lqStats && step == 0) {
                const IntSpan fp0 = wallBucket.pointsOf(fi);
                std::vector<Vec3> loop0;
                loop0.reserve(static_cast<std::size_t>(fp0.size()));
                for (int p : fp0) loop0.push_back(out.points[static_cast<std::size_t>(p)]);
                lqSeedNonPlanarity[fi] = loopNonPlanarity(loop0);
            }
            if (reverted[static_cast<std::size_t>(fi)]) {
                newWallBucket.appendFrom(wallBucket, fi); // untouched revert
                newFaceOrigin.push_back(faceOrigin[static_cast<std::size_t>(fi)]);
                newFaceGateStack.push_back(faceGateStack[static_cast<std::size_t>(fi)]);
                continue;
            }
            // A merged face extends the prism above it (see mergeInto): no new
            // cell, and its top face -- the interface -- is not emitted.
            const bool merged = mergeInto[static_cast<std::size_t>(fi)] >= 0;
            const int prismCell = merged ? mergeInto[static_cast<std::size_t>(fi)] : nextCellIndex++;
            if (!merged) ++stats.perStepPrismCells[static_cast<std::size_t>(step)];
            topFaceToPrism[static_cast<std::size_t>(fi)] = prismCell;
            const IntSpan topPts = wallBucket.pointsOf(fi);
            for (int p : topPts) {
                if (static_cast<std::size_t>(p) >= pointSpent.size()) pointSpent.resize(out.points.size(), 0);
                pointSpent[static_cast<std::size_t>(p)] = 1;
            }
            const int coreCell = wallBucket.owner[static_cast<std::size_t>(fi)];
            if (!merged) {
                if (static_cast<std::size_t>(prismCell) >= cellReached.size()) cellReached.resize(static_cast<std::size_t>(prismCell) + 1, 1);
                cellReached[static_cast<std::size_t>(prismCell)] =
                    static_cast<std::size_t>(coreCell) < cellReached.size() ? cellReached[static_cast<std::size_t>(coreCell)] : 1;
                // Top face becomes internal: owner (core) < neighbour
                // (prism) always holds since prism indices are appended
                // after every pre-existing cell.
                internalOut.appendFrom(wallBucket, fi, coreCell, prismCell, -1);
            }

            std::vector<int> bottomPts;
            bottomPts.reserve(static_cast<std::size_t>(topPts.size()));
            for (int p : topPts) {
                bottomPts.push_back(bottomPointIdx[static_cast<std::size_t>(frontIdx[p])]);
            }
            newWallBucket.append(bottomPts, prismCell, -1, pOrd);
            newFaceOrigin.push_back(faceOrigin[static_cast<std::size_t>(fi)]);
            newFaceGateStack.push_back(faceGateStack[static_cast<std::size_t>(fi)]);
            if (merged) {
                mergedPrisms.insert(prismCell);
            } else {
                prismTopOf[prismCell] = std::vector<int>(topPts.begin(), topPts.end());
                prismBottomOf[prismCell] = bottomPts;
                prismParentOf[prismCell] = coreCell;
            }

            // Evaluate this newly-created cell's own quality.
            // ALWAYS on (this is the gate's predicate; cost is
            // per newly-created prism only). Origin id used to key the
            // flagged-stack set is the ORIGIN this face's ancestry
            // chains back to (minted at step 0), not this step's local
            // fi -- a stack is condemned as a whole, across steps.
            {
                // Evaluated once, in the validation loop above (a face
                // reaching emission has passed it); reused here for the
                // next step's owner centroid and the safety-net report.
                const LayerCellQuality& q = faceQuality[static_cast<std::size_t>(fi)];
                lqPrismCentroid[prismCell] = q.centroid;
                if (lqStats) {
                    lqWorstPerCell.push_back(q.worstTetNorm);
                    lqBadPyramidFaces += q.nBadPyramidFaces;
                    lqBadTetFaces += q.nBadTetFaces;
                }
                if (q.nBadPyramidFaces > 0 || q.nBadTetFaces > 0 || q.nHighSkewFaces > 0 || q.nHighNonOrthFaces > 0) {
                    ++lqFlaggedCells;
                    const int origin = faceGateStack[static_cast<std::size_t>(fi)];
                    if (origin < 0) {
                        // No stack to condemn (not reachable while every
                        // wall face -- seams included, see faceGateStack
                        // -- carries a gate identity; kept as a loud
                        // disclosure rather than an assert).
                        ++report.orphanBadCells;
                    } else if (lqFlaggedStacks.insert({pOrd, origin}).second) {
                        report.badStacks.insert({pOrd, origin});
                        if (lqStats) {
                            const auto it = lqSeedNonPlanarity.find(origin);
                            const double nonPlanarity = it != lqSeedNonPlanarity.end() ? it->second : -1.0;
                            lqNewFlaggedThisStep.emplace_back(origin, nonPlanarity);
                        }
                        if (gateDump) {
                            const int coreLevelNow = levelOfCell(coreCell);
                            report.condemnedSites.push_back(
                                {pOrd, origin, q.centroid.x, q.centroid.y, q.centroid.z, coreLevelNow});
                        }
                    }
                }
            }

            const int coreLevel = levelOfCell(coreCell);
            if (!merged) {
                newCellLevels.push_back(coreLevel);
                // Tag this prism with the stack it belongs to and the
                // march step that made it (parallel to newCellLevels, so the same
                // append order and the same final compaction apply).
                newCellStackId.push_back(faceOrigin[static_cast<std::size_t>(fi)]);
                newCellLayerIndex.push_back(step);
            }
            for (int p : topPts) {
                const std::size_t di = static_cast<std::size_t>(frontIdx[p]);
                bottomPointLevelAcc[di] = std::max(bottomPointLevelAcc[di], coreLevel);
            }

            // Residual metric: distance-to-STL at this wall face's
            // (post-landing) centroid, area-weighted. LAST layer step
            // only -- intermediate fronts sit at layer interfaces by
            // design, their distance to the STL is not a residual.
            // (Faces dropped on an earlier step keep their reverted
            // position and are NOT counted -- droppedFaces reports
            // them separately.)
            if (lastStep) {
                std::vector<Vec3> botLoop;
                botLoop.reserve(bottomPts.size());
                for (int p : bottomPts) botLoop.push_back(out.points[static_cast<std::size_t>(p)]);
                const Vec3 cen = loopCentroid(botLoop);
                const double area = loopArea(botLoop);
                ClosestHit h = closestPointOnSoup(bins, cen);
                const double res = std::sqrt(h.distSq);
                residualAreaSum += area;
                residualWeightedSum += area * res;
                residualMax = std::max(residualMax, res);
            }

            // Side faces, shared between adjacent kept stacks via the
            // shared edge key (built once per edge, reused).
            const int n = topPts.size();
            for (int i = 0; i < n; ++i) {
                const int a = topPts[i];
                const int b = topPts[(i + 1) % n];
                const EdgeKey key = makeEdgeKey(a, b);
                // Find the OTHER top face sharing this edge (if any,
                // and if it is also kept) to decide whether the side
                // face is interior or a boundary.
                int otherFi = -1;
                auto it = edgeToFaces.find(key);
                if (it != edgeToFaces.end()) {
                    for (int cand : it->second) {
                        if (cand != fi) { otherFi = cand; break; }
                    }
                }
                const bool otherKept = otherFi != -1 && !reverted[static_cast<std::size_t>(otherFi)];
                if (otherKept && otherFi < fi) {
                    // Already emitted by the other (lower-index) prism
                    // when IT processed this edge -- just needs its
                    // neighbour set, done below in a second pass.
                    continue;
                }
                // Outward-from-THIS-prism winding is fixed by topology:
                // the top loop (outward-from-core, i.e. pointing INTO the
                // prism) traverses a->b, so the prism's outward side quad
                // is {a, b, bot(b), bot(a)} -- always. Never re-wound
                // geometrically (the pre-fix flip produced open cells,
                // see prismValid).
                const std::vector<int> quad{a, b, bottomPointIdx[static_cast<std::size_t>(frontIdx[b])],
                                             bottomPointIdx[static_cast<std::size_t>(frontIdx[a])]};
                if (otherKept) {
                    // otherFi > fi (we skipped the < case above): this
                    // face becomes internal once the other prism is
                    // created; record it as a pending internal face by
                    // pushing now with neighbour resolved in the second
                    // pass below via a lookup keyed by edge.
                    internalOut.append(quad, prismCell, -1, -1); // neighbour patched in the pass below
                    pendingSideByEdge[key].push_back(internalOut.size() - 1);
                } else {
                    const std::vector<int> planesA = domainPlanesOf(out.points[static_cast<std::size_t>(a)], cfg);
                    const std::vector<int> planesB = domainPlanesOf(out.points[static_cast<std::size_t>(b)], cfg);
                    int commonPlane = -1;
                    for (int pl : planesA) {
                        if (std::find(planesB.begin(), planesB.end(), pl) != planesB.end()) { commonPlane = pl; break; }
                    }
                    if (commonPlane >= 0) {
                        const int domPatchOrd = patchNameToOrdinal.at(domainPatchNameOf(commonPlane, cfg));
                        boundaryByPatch[static_cast<std::size_t>(domPatchOrd)].append(quad, prismCell, -1,
                                                                                      domPatchOrd);
                    } else {
                        // Seam against a dropped/un-emitted neighbour (or
                        // any other boundary case): closes the region as
                        // a wall face of the SAME solid, per the design's
                        // drop-locally seam-closure rule.
                        newWallBucket.append(quad, prismCell, -1, pOrd);
                        newFaceOrigin.push_back(-1); // seam face: no original wall-face ancestor
                        newFaceGateStack.push_back(faceGateStack[static_cast<std::size_t>(fi)]);
                    }
                }
            }
        }
        // Second pass: patch the pending shared-side-face neighbours.
        // These were appended to internalOut in creation order right
        // after each prism's own top-internal face; find them again by
        // owner+points signature (small counts per case, fine to scan).
        for (int fi = 0; fi < nTop; ++fi) {
            if (reverted[static_cast<std::size_t>(fi)]) continue;
            const IntSpan topPts = wallBucket.pointsOf(fi);
            const int n = topPts.size();
            for (int i = 0; i < n; ++i) {
                const int a = topPts[i];
                const int b = topPts[(i + 1) % n];
                auto it = edgeToFaces.find(makeEdgeKey(a, b));
                if (it == edgeToFaces.end()) continue;
                int otherFi = -1;
                for (int cand : it->second) {
                    if (cand != fi) { otherFi = cand; break; }
                }
                if (otherFi == -1 || reverted[static_cast<std::size_t>(otherFi)] || otherFi >= fi) continue;
                // fi is the higher-index prism sharing this edge with
                // otherFi (< fi), which already emitted the shared side
                // face with neighbour left at -1 -- set it now.
                // Candidates come from the per-edge pending index (see
                // the index's comment above); the predicate does not
                // depend on the index.
                auto pendIt = pendingSideByEdge.find(makeEdgeKey(a, b));
                if (pendIt == pendingSideByEdge.end()) continue;
                for (int cand : pendIt->second) {
                    if (internalOut.owner[static_cast<std::size_t>(cand)] ==
                            topFaceToPrism[static_cast<std::size_t>(otherFi)] &&
                        internalOut.neighbour[static_cast<std::size_t>(cand)] == -1 &&
                        internalOut.pointCount(cand) == 4) {
                        std::set<int> want{a, b, bottomPointIdx[static_cast<std::size_t>(frontIdx[a])],
                                            bottomPointIdx[static_cast<std::size_t>(frontIdx[b])]};
                        const IntSpan cp = internalOut.pointsOf(cand);
                        std::set<int> have(cp.begin(), cp.end());
                        if (want == have) {
                            internalOut.neighbour[static_cast<std::size_t>(cand)] =
                                topFaceToPrism[static_cast<std::size_t>(fi)];
                            // TWO-SIDED test, now that both prisms exist. The
                            // per-cell evaluation above scores a side quad
                            // against the neighbour's CANDIDATE prism, and its
                            // tet test against its own centre only; checkMesh
                            // scores an INTERNAL face against BOTH final centres
                            // (findSharedBasePoint). MEASURED (win_wfp_apron with
                            // the per-vertex half-grid band): 5 internal side
                            // quads between adjacent triangular prisms shipped
                            // with "Error in face tets" through that gap. Both
                            // stacks are condemned -- the face belongs to both.
                            std::vector<Vec3> quadLoop;
                            quadLoop.reserve(4);
                            for (int p : cp) quadLoop.push_back(out.points[static_cast<std::size_t>(p)]);
                            const Vec3& ownC = faceQuality[static_cast<std::size_t>(otherFi)].centroid;
                            const Vec3& nbrC = faceQuality[static_cast<std::size_t>(fi)].centroid;
                            if (!faceHasUsableBasePoint(ownC, &nbrC, quadLoop, kMinTetQuality) ||
                                faceNonOrthDegOf(quadLoop, ownC, &nbrC) > kNonOrthMaxDeg ||
                                faceSkewnessOf(quadLoop, ownC, &nbrC) > kSkewMax) {
                                for (int f2 : {otherFi, fi}) {
                                    const int origin2 = faceGateStack[static_cast<std::size_t>(f2)];
                                    if (origin2 < 0) {
                                        ++report.orphanBadCells;
                                    } else if (lqFlaggedStacks.insert({pOrd, origin2}).second) {
                                        report.badStacks.insert({pOrd, origin2});
                                    }
                                }
                            }
                            break;
                        }
                    }
                }
            }
        }

        for (std::size_t i = 0; i < nFront; ++i) {
            pointLevel[static_cast<std::size_t>(bottomPointIdx[i])] =
                bottomPointLevelAcc[i] >= 0 ? bottomPointLevelAcc[i] : 0;
        }

        // Per-march-step layer-quality report (opt-in,
        // NINJA_LAYER_QUALITY_STATS). Mimics Cutter.cpp's "wall-facet
        // stats" block style/delimiters.
        if (lqStats) {
            std::cout << "--- layer-quality stats (patch " << spec.wallPatchName << ", layer " << layerJ << ") ---\n";
            std::cout << "new cells this layer                = "
                      << stats.perStepPrismCells[static_cast<std::size_t>(step)] << "\n";
            std::cout << "flagged cells (pyramid and/or tet)   = " << lqFlaggedCells << "\n";
            std::cout << "flagged faces: pyramid orientation   = " << lqBadPyramidFaces << "\n";
            std::cout << "flagged faces: tet quality           = " << lqBadTetFaces << "\n";
            std::cout << "worst-cell tet-quality histogram (normalized by layer-step^3):\n";
            {
                // Bucket edges chosen around the calibrated 1e-9 flag
                // threshold, log-scale outward in both directions.
                static const double edges[] = {-1e-3, -1e-6, -1e-9, 0.0, 1e-9, 1e-6, 1e-3};
                const int nBuckets = 8;
                std::vector<int> hist(static_cast<std::size_t>(nBuckets), 0);
                for (double w : lqWorstPerCell) {
                    int b = 0;
                    while (b < 6 && w >= edges[b]) ++b;
                    ++hist[static_cast<std::size_t>(b)];
                }
                const char* labels[8] = {"< -1e-3", "[-1e-3,-1e-6)", "[-1e-6,-1e-9)", "[-1e-9,0)",
                                          "[0,1e-9)", "[1e-9,1e-6)", "[1e-6,1e-3)", ">= 1e-3"};
                for (int b = 0; b < nBuckets; ++b) {
                    std::cout << "  " << std::left << std::setw(16) << labels[b] << " " << std::right << std::setw(8)
                              << hist[static_cast<std::size_t>(b)] << "\n";
                }
            }
            if (!lqNewFlaggedThisStep.empty()) {
                std::cout << "newly flagged stacks this layer (origin id, seed non-planarity):\n";
                for (const auto& [origin, nonPlanarity] : lqNewFlaggedThisStep) {
                    std::cout << "  origin=" << std::left << std::setw(8) << origin << std::right
                              << " nonPlanarity=" << nonPlanarity << "\n";
                }
            }
            std::cout << "--- end layer-quality stats ---\n";
        }

        // Progress line, flushed immediately -- the march is the long
        // stage on big cases and must show signs of life per step.
        std::cout << "layers " << spec.wallPatchName << ": step " << (step + 1) << "/" << nSteps << ", "
                  << stats.perStepPrismCells[static_cast<std::size_t>(step)] << " prism cells" << std::endl;

        wallBucket = std::move(newWallBucket);
        faceOrigin = std::move(newFaceOrigin);
        faceGateStack = std::move(newFaceGateStack);
        clk.lap("emission");
        } // step loop (layer-by-layer march)

        // Propagate this STL's frozen/dropped origin sets onto the
        // actual final mesh point ids. `frozenOrigins` is
        // POINT-origin-keyed (`originOfPoint`), so it maps onto
        // `pointOrigin` directly; `droppedOrigins` is FACE-origin-keyed
        // (`faceOrigin[fi]`, a DIFFERENT id space -- conflating the two
        // leaves `dropped` almost entirely unset even when many faces
        // dropped). Dropped points are therefore derived from the
        // CURRENT (post-step-loop) `wallBucket` + `faceOrigin` directly:
        // every point of every face whose own origin is in
        // `droppedOrigins`. Grows the flag arrays to the current point
        // count -- `out.points` only ever grows (bottom points appended
        // per step).
        if (pointFrozenFlag.size() < out.points.size()) pointFrozenFlag.resize(out.points.size(), 0);
        if (pointDroppedFlag.size() < out.points.size()) pointDroppedFlag.resize(out.points.size(), 0);
        for (const auto& [pt, origin] : pointOrigin) {
            if (pt < 0 || static_cast<std::size_t>(pt) >= pointFrozenFlag.size()) continue;
            if (frozenOrigins.count(origin)) pointFrozenFlag[static_cast<std::size_t>(pt)] = 1;
        }
        for (int fi = 0; fi < wallBucket.size(); ++fi) {
            if (!droppedOrigins.count(faceOrigin[static_cast<std::size_t>(fi)])) continue;
            for (int p : wallBucket.pointsOf(fi)) {
                if (p < 0 || static_cast<std::size_t>(p) >= pointDroppedFlag.size()) continue;
                pointDroppedFlag[static_cast<std::size_t>(p)] = 1;
            }
        }
    }
    (void)layeredPatchOrdinals;

    report.totalWallFaces0 = lqTotalWallFaces0;

    // Per-pass summary (opt-in,
    // NINJA_LAYER_QUALITY_STATS). "flagged" here is what THIS pass newly
    // condemns; the cumulative REMOVED count is reported by applyLayers.
    if (lqStats) {
        const double frac = lqTotalWallFaces0 > 0
                                 ? static_cast<double>(lqFlaggedStacks.size()) / static_cast<double>(lqTotalWallFaces0)
                                 : 0.0;
        std::cout << "--- layer-quality stats (pass summary) ---\n";
        std::cout << "total wall seed faces               = " << lqTotalWallFaces0 << "\n";
        std::cout << "stacks already removed on entry     = " << gateRemoved.size() << "\n";
        std::cout << "newly flagged stacks                = " << lqFlaggedStacks.size() << "\n";
        std::cout << "newly flagged fraction              = " << frac << "\n";
        std::cout << "flagged cells with no removable stack = " << report.orphanBadCells << "\n";
        std::cout << "--- end layer-quality stats ---\n";
    }

    clk.lap("spec tail");
    // --- Reassemble: enforce owner < neighbour on every internal face
    // (swap + reverse otherwise), then sort by (owner, neighbour).
    for (int fi = 0; fi < internalOut.size(); ++fi) {
        int& fo = internalOut.owner[static_cast<std::size_t>(fi)];
        int& fn = internalOut.neighbour[static_cast<std::size_t>(fi)];
        if (fn != -1 && fo > fn) {
            std::swap(fo, fn);
            std::reverse(internalOut.mutablePointsOf(fi),
                         internalOut.mutablePointsOf(fi) + internalOut.pointCount(fi));
        }
    }
    internalOut.stableSortByOwnerNeighbour();

    GeneratedMesh finalMesh;
    finalMesh.points = out.points;
    finalMesh.nInternalFaces = internalOut.size();
    finalMesh.faces = std::move(internalOut);
    for (std::size_t p = 0; p < out.patches.size(); ++p) {
        PatchInfo info;
        info.name = out.patches[p].name;
        info.type = out.patches[p].type;
        info.startFace = finalMesh.faces.size();
        info.nFaces = boundaryByPatch[p].size();
        finalMesh.faces.appendAll(boundaryByPatch[p]);
        boundaryByPatch[p].clear();
        finalMesh.patches.push_back(info);
    }

    const int nCellsTotal = nextCellIndex;
    buildCellFaces(finalMesh, nCellsTotal);

    cellLevel.insert(cellLevel.end(), newCellLevels.begin(), newCellLevels.end());
    {
        // -1 for every pre-existing (non-prism) cell, then the per-prism tags
        // in the same append order.
        std::vector<int>& sid = result.cellStackId;
        std::vector<int>& lix = result.cellLayerIndex;
        sid.assign(cellLevelIn.size(), -1);
        lix.assign(cellLevelIn.size(), -1);
        sid.insert(sid.end(), newCellStackId.begin(), newCellStackId.end());
        lix.insert(lix.end(), newCellLayerIndex.begin(), newCellLayerIndex.end());
    }

    clk.lap("reassemble");
    // --- Final point compaction: some bottom points are created for
    // EVERY front point (needed to evaluate prism validity for dropped
    // candidates too) but a point whose only incident faces all got
    // dropped is never referenced by the final face set -- mirrors
    // Cutter.cpp's own compaction pass (ascending-original-index order,
    // exact-position dedup, same "bit-identical keying only" rule).
    {
        std::set<int> used(finalMesh.faces.points.begin(), finalMesh.faces.points.end());
        std::unordered_map<int, int> remap;
        std::vector<Vec3> compactPoints;
        std::vector<int> compactLevel;
        std::map<std::tuple<double, double, double>, int> coordToNew;
        compactPoints.reserve(used.size());
        // OR the per-node diagnostic flags through the
        // exact same dedup remap (a coordinate collapse should carry
        // either contributor's flag).
        std::vector<char> compactFrozen, compactDropped;
        auto flagAt = [](const std::vector<char>& flags, int idx) -> char {
            return static_cast<std::size_t>(idx) < flags.size() ? flags[static_cast<std::size_t>(idx)] : 0;
        };
        for (int oldIdx : used) {
            const Vec3& v = finalMesh.points[static_cast<std::size_t>(oldIdx)];
            const auto key = std::make_tuple(v.x, v.y, v.z);
            auto it = coordToNew.find(key);
            if (it != coordToNew.end()) {
                remap[oldIdx] = it->second;
                compactLevel[static_cast<std::size_t>(it->second)] =
                    std::max(compactLevel[static_cast<std::size_t>(it->second)], pointLevel[static_cast<std::size_t>(oldIdx)]);
                compactFrozen[static_cast<std::size_t>(it->second)] =
                    compactFrozen[static_cast<std::size_t>(it->second)] || flagAt(pointFrozenFlag, oldIdx);
                compactDropped[static_cast<std::size_t>(it->second)] =
                    compactDropped[static_cast<std::size_t>(it->second)] || flagAt(pointDroppedFlag, oldIdx);
                continue;
            }
            const int newIdx = static_cast<int>(compactPoints.size());
            coordToNew.emplace(key, newIdx);
            remap[oldIdx] = newIdx;
            compactPoints.push_back(v);
            compactLevel.push_back(pointLevel[static_cast<std::size_t>(oldIdx)]);
            compactFrozen.push_back(flagAt(pointFrozenFlag, oldIdx));
            compactDropped.push_back(flagAt(pointDroppedFlag, oldIdx));
        }
        // Remap + collapse-consecutive-duplicates, streaming into a
        // fresh CSR store (the collapse can SHORTEN a loop, so the face
        // list is rebuilt rather than edited in place -- same face
        // order, same per-face result as the former in-place loop).
        // Landing can move TWO front points onto the SAME closest
        // point (measured, offset_rotcube: adjacent front points
        // flanking the cube's corner crease at equal height both
        // snap to the identical crease-line point). The exact-
        // position dedup above then correctly fuses them into one
        // point -- but the face loop is left with the same index
        // twice in a row (a quad that is geometrically a triangle).
        // Collapse consecutive duplicates (wrap-around included);
        // loops that would fall below 3 points are left untouched
        // (disclosed fallback -- integrity would flag them, none
        // observed).
        FaceStore rebuilt;
        rebuilt.points.reserve(finalMesh.faces.points.size());
        rebuilt.offsets.reserve(finalMesh.faces.offsets.size());
        rebuilt.owner.reserve(finalMesh.faces.owner.size());
        rebuilt.neighbour.reserve(finalMesh.faces.neighbour.size());
        rebuilt.patchId.reserve(finalMesh.faces.patchId.size());
        std::vector<int> loopScratch, collapsed;
        for (int fi = 0; fi < finalMesh.faces.size(); ++fi) {
            const IntSpan fp = finalMesh.faces.pointsOf(fi);
            loopScratch.clear();
            for (int p : fp) loopScratch.push_back(remap.at(p));
            collapsed.clear();
            for (std::size_t i = 0; i < loopScratch.size(); ++i) {
                if (collapsed.empty() || loopScratch[i] != collapsed.back()) {
                    collapsed.push_back(loopScratch[i]);
                }
            }
            while (collapsed.size() > 1 && collapsed.front() == collapsed.back()) {
                collapsed.pop_back();
            }
            const std::vector<int>& emit = collapsed.size() >= 3 ? collapsed : loopScratch;
            rebuilt.append(emit, finalMesh.faces.owner[static_cast<std::size_t>(fi)],
                           finalMesh.faces.neighbour[static_cast<std::size_t>(fi)],
                           finalMesh.faces.patchId[static_cast<std::size_t>(fi)]);
        }
        finalMesh.faces = std::move(rebuilt);
        finalMesh.points = std::move(compactPoints);
        pointLevel = std::move(compactLevel);
        pointFrozenFlag = std::move(compactFrozen);
        pointDroppedFlag = std::move(compactDropped);
    }

    stats.residualMean = residualAreaSum > 0.0 ? residualWeightedSum / residualAreaSum : 0.0;
    stats.residualMax = residualMax;

    clk.lap("compaction");
    clk.print();

    result.mesh = std::move(finalMesh);
    result.cellLevel = std::move(cellLevel);
    result.pointLevel = std::move(pointLevel);
    result.stats = stats;
    report.sealedArea = stats.sealDroppedArea;
    report.thinArea = stats.thinDroppedArea;
    result.pointFrozen = std::move(pointFrozenFlag);
    result.pointDropped = std::move(pointDroppedFlag);
    return result;
}

} // namespace

long reattributeSteepWallFaces(GeneratedMesh& mesh, const std::vector<LayerStlSpec>& specs,
                               const std::vector<TriangleAabbBins>& perStlBins) {
    if (specs.size() < 2) return 0;
    std::unordered_map<std::string, int> ordOf;
    for (std::size_t p = 0; p < mesh.patches.size(); ++p) ordOf[mesh.patches[p].name] = static_cast<int>(p);
    std::vector<int> specOfPatch(mesh.patches.size(), -1);
    for (std::size_t k = 0; k < specs.size(); ++k) specOfPatch[static_cast<std::size_t>(ordOf.at(specs[k].wallPatchName))] = static_cast<int>(k);
    // Cosine of the face normal with the direction to STL k's closest point,
    // and the face centre's offset error |d_k - t_k| / t_k.
    auto facing = [&](const Vec3& n, const Vec3& c, std::size_t k, double& relErr) {
        const ClosestHit h = closestPointOnSoup(perStlBins[static_cast<std::size_t>(specs[k].stlIndex)], c);
        const double tk = specs[k].localThickness ? specs[k].localThickness(c) : specs[k].thickness;
        relErr = std::fabs(std::sqrt(h.distSq) - tk) / tk;
        return dot(n, normalizeOrZero(h.point - c));
    };
    const int nB = mesh.nFaces() - mesh.nInternalFaces;
    std::vector<int> newPatch(static_cast<std::size_t>(nB), -1);
    for (std::size_t p = 0; p < mesh.patches.size(); ++p)
        for (int f = mesh.patches[p].startFace; f < mesh.patches[p].startFace + mesh.patches[p].nFaces; ++f)
            newPatch[static_cast<std::size_t>(f - mesh.nInternalFaces)] = static_cast<int>(p);
    long moved = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256) reduction(+ : moved)
#endif
    for (int b = 0; b < nB; ++b) {
        const int f = mesh.nInternalFaces + b;
        const int p = newPatch[static_cast<std::size_t>(b)];
        if (p < 0) continue;
        const int own = specOfPatch[static_cast<std::size_t>(p)];
        if (own < 0) continue;
        std::vector<Vec3> loop;
        for (int q : mesh.faces.pointsOf(f)) loop.push_back(mesh.points[static_cast<std::size_t>(q)]);
        const Vec3 n = normalizeOrZero(newellNormal(loop));
        const Vec3 c = loopCentroid(loop);
        double err;
        const double cosOwn = facing(n, c, static_cast<std::size_t>(own), err);
        if (cosOwn >= 0.5) continue;
        double bestCos = std::max(cosOwn, 0.8);
        int best = -1;
        for (std::size_t k = 0; k < specs.size(); ++k) {
            if (static_cast<int>(k) == own) continue;
            const double ck = facing(n, c, k, err);
            if (err <= 0.5 && ck > bestCos) {
                bestCos = ck;
                best = ordOf.at(specs[k].wallPatchName);
            }
        }
        if (best >= 0) {
            newPatch[static_cast<std::size_t>(b)] = best;
            ++moved;
        }
    }
    if (moved == 0) return 0;
    // Regroup the boundary faces by patch, each patch keeping its face order.
    FaceStore ns;
    for (int f = 0; f < mesh.nInternalFaces; ++f) ns.appendFrom(mesh.faces, f);
    for (std::size_t p = 0; p < mesh.patches.size(); ++p) {
        mesh.patches[p].startFace = ns.size();
        for (int b = 0; b < nB; ++b) {
            if (newPatch[static_cast<std::size_t>(b)] != static_cast<int>(p)) continue;
            const int f = mesh.nInternalFaces + b;
            const IntSpan fp = mesh.faces.pointsOf(f);
            ns.append(fp.b, fp.size(), mesh.faces.owner[static_cast<std::size_t>(f)], -1, static_cast<int>(p));
        }
        mesh.patches[p].nFaces = ns.size() - mesh.patches[p].startFace;
    }
    mesh.faces = std::move(ns);
    buildCellFaces(mesh, mesh.nCells());
    return moved;
}

LayersResult applyLayers(const GeneratedMesh& cutMeshIn, const std::vector<int>& cellLevelIn,
                          const std::vector<int>& pointLevelIn, const MeshConfig& cfg,
                          const std::vector<LayerStlSpec>& specs,
                          const std::vector<std::vector<Triangle>>& perStlTris,
                          const std::vector<TriangleAabbBins>& perStlBins,
                          const std::vector<ForceDropSphere>& forceDropSpheres, const Vec3& locationInMesh) {
    // --- The offset cut's own disclosure, measured on the mesh the
    // march is HANDED, before it touches anything.
    //
    // MEASURED (bm_layers_descargador): the layer stage's input has
    // 354,800 cells / 1700 m3 of fluid in 17 components that cannot
    // reach locationInMesh, against 108,259 cells / 1018 m3 (ONE
    // geometric feature, the rock void the case dict documents) in the
    // unlayered sibling. The extra ~680 m3 is not the march's and not
    // the smoothing's: cutting the core at the offset level set d = t
    // removes a shell of thickness t from EVERY wall, so every fluid
    // passage narrower than 2t is gone before a single prism exists,
    // and each pocket behind such a passage is discarded by the final
    // flood fill. The march cannot give it back -- there are no cells
    // in the passage to march from -- so the only honest thing this
    // stage can do is say so, loudly, with the numbers, every run.
    // (The remedy is a smaller total thickness, or local refinement:
    // 2t must stay under the narrowest passage that carries flow.)
    {
        const std::vector<char> reached = reachableCellsFromLocation(cutMeshIn, locationInMesh);
        const int nCells = cutMeshIn.nCells();
        // Per-cell volume by the divergence formula over its faces,
        // each face fanned about its own point average -- the same
        // decomposition loopArea/loopCentroid use elsewhere here.
        std::vector<double> cellVol(static_cast<std::size_t>(nCells), 0.0);
        const int nF = cutMeshIn.nFaces();
        for (int f = 0; f < nF; ++f) {
            const IntSpan fp = cutMeshIn.faces.pointsOf(f);
            const int n = fp.size();
            if (n < 3) continue;
            std::vector<Vec3> loop;
            loop.reserve(static_cast<std::size_t>(n));
            for (int p : fp) loop.push_back(cutMeshIn.points[static_cast<std::size_t>(p)]);
            const Vec3 area = newellNormal(loop) * 0.5;
            const Vec3 c = loopCentroid(loop);
            const double contrib = dot(area, c) / 3.0;
            const int own = cutMeshIn.faces.owner[static_cast<std::size_t>(f)];
            const int nei = cutMeshIn.faces.neighbour[static_cast<std::size_t>(f)];
            if (own >= 0) cellVol[static_cast<std::size_t>(own)] += contrib;
            if (nei >= 0) cellVol[static_cast<std::size_t>(nei)] -= contrib;
        }
        // Connected components over the UNREACHED cells only.
        std::vector<int> comp(static_cast<std::size_t>(nCells), -1);
        int nComp = 0;
        double lostVol = 0.0;
        long lostCells = 0;
        struct Lost {
            long cells = 0;
            double vol = 0.0;
            Vec3 lo{0, 0, 0}, hi{0, 0, 0};
        };
        std::vector<Lost> lost;
        std::vector<int> stack;
        for (int c0 = 0; c0 < nCells; ++c0) {
            if (reached[static_cast<std::size_t>(c0)] || comp[static_cast<std::size_t>(c0)] >= 0) continue;
            const int id = nComp++;
            Lost rec;
            bool first = true;
            comp[static_cast<std::size_t>(c0)] = id;
            stack.clear();
            stack.push_back(c0);
            while (!stack.empty()) {
                const int c = stack.back();
                stack.pop_back();
                ++rec.cells;
                rec.vol += cellVol[static_cast<std::size_t>(c)];
                for (int f : cutMeshIn.cellFacesOf(c)) {
                    for (int p : cutMeshIn.faces.pointsOf(f)) {
                        const Vec3& v = cutMeshIn.points[static_cast<std::size_t>(p)];
                        if (first) {
                            rec.lo = v;
                            rec.hi = v;
                            first = false;
                        } else {
                            rec.lo = Vec3{std::min(rec.lo.x, v.x), std::min(rec.lo.y, v.y), std::min(rec.lo.z, v.z)};
                            rec.hi = Vec3{std::max(rec.hi.x, v.x), std::max(rec.hi.y, v.y), std::max(rec.hi.z, v.z)};
                        }
                    }
                    const int own = cutMeshIn.faces.owner[static_cast<std::size_t>(f)];
                    const int nei = cutMeshIn.faces.neighbour[static_cast<std::size_t>(f)];
                    const int other = own == c ? nei : own;
                    if (other >= 0 && !reached[static_cast<std::size_t>(other)] &&
                        comp[static_cast<std::size_t>(other)] < 0) {
                        comp[static_cast<std::size_t>(other)] = id;
                        stack.push_back(other);
                    }
                }
            }
            lostCells += rec.cells;
            lostVol += rec.vol;
            lost.push_back(rec);
        }
        if (nComp > 0) {
            std::sort(lost.begin(), lost.end(), [](const Lost& a, const Lost& b) { return a.vol > b.vol; });
            std::cout << "offsetCutDisconnectedComponents = " << nComp << "\n"
                      << "offsetCutDisconnectedCells = " << lostCells << "\n"
                      << "offsetCutDisconnectedVolume = " << lostVol << "\n";
            const std::size_t nShow = std::min<std::size_t>(lost.size(), 10);
            for (std::size_t i = 0; i < nShow; ++i) {
                std::cout << "offsetCutDisconnected[" << i << "] cells = " << lost[i].cells
                          << " volume = " << lost[i].vol << " bbox = (" << lost[i].lo.x << " " << lost[i].lo.y
                          << " " << lost[i].lo.z << ")..(" << lost[i].hi.x << " " << lost[i].hi.y << " "
                          << lost[i].hi.z << ")\n";
            }
            std::cerr << "warning: layers: the offset cut (d = t) left " << lostVol << " m3 of fluid in " << nComp
                      << " component(s) unreachable from locationInMesh; the final flood fill will DISCARD it."
                      << " Every passage narrower than 2*t is closed before the march runs -- reduce the total"
                      << " layer thickness (or refine locally) if that fluid carries flow\n";
        }
    }

    // The gate's outer loop. Design note (why a
    // re-march instead of an in-march retraction): a stack is only known
    // to be bad AFTER at least one of its prisms exists, and the march
    // has no "un-create a cell" primitive -- unwinding one would mean
    // rolling back appended points/cells/faces AND the front the
    // neighbouring stacks were already marched against, i.e. exactly the
    // post-hoc surgery the design forbids. Pre-loading the condemned set
    // and re-marching instead reuses the drop path verbatim, so a
    // removed stack is INDISTINGUISHABLE from a face the march never
    // extruded. The condemned set only grows, one pass strictly
    // dominates the last, and nothing depends on iteration order
    // (`std::set`-keyed by (patch ordinal, origin id)), so the result is
    // deterministic and the loop terminates.
    //
    // Passes are capped: removing a stack perturbs its neighbours'
    // march, so a pass can uncover NEW bad stacks (the "removal
    // closure"). In practice this settles in very few passes; the cap only
    // bounds pathological cases, and a residual is disclosed loudly
    // rather than silently accepted.
    constexpr int kMaxGatePasses = 6;
    std::set<std::pair<int, int>> removed;
    std::size_t initiallyBad = 0;
    LayersResult result;
    LayerGateReport report;
    int passesRun = 0;
    // Union of condemned-stack sites over every pass,
    // deduplicated by (patch ordinal, origin) so a stack rediscovered on
    // a later pass (e.g. a terrace-closure neighbour) is only recorded
    // once, at its first flagging.
    const bool gateDumpOuter = std::getenv("NINJA_LAYER_GATE_DUMP") != nullptr;
    // Seal-guard disclosure, accumulated over every pass. A stack the
    // guard condemns is gone from the NEXT pass's landing, so the final
    // pass's own sealed-region instrument reads zero -- which is the
    // right mesh but the wrong report. These carry the fact across.
    std::set<std::pair<int, int>> sealedAll;
    double sealedAreaAll = 0.0;
    long orphanSealedAll = 0;
    std::set<std::pair<int, int>> thinAll;
    double thinAreaAll = 0.0;
    long orphanThinAll = 0;
    std::vector<CondemnedSite> allCondemnedSites;
    std::set<std::pair<int, int>> dumpedKeys;
    for (int pass = 0; pass < kMaxGatePasses; ++pass) {
        report = LayerGateReport{};
        result = applyLayersPass(cutMeshIn, cellLevelIn, pointLevelIn, cfg, specs, perStlTris, perStlBins,
                                 forceDropSpheres, removed, report, locationInMesh);
        ++passesRun;
        std::cout << "layers: gate pass " << passesRun << ", " << report.badStacks.size()
                  << " newly condemned stacks" << std::endl;
        if (pass == 0) initiallyBad = report.badStacks.size();
        if (gateDumpOuter) {
            for (const CondemnedSite& s : report.condemnedSites) {
                if (dumpedKeys.insert({s.patchOrdinal, s.origin}).second) {
                    allCondemnedSites.push_back(s);
                }
            }
        }
        if (!report.sealedStacks.empty() || report.orphanSealedFaces > 0) {
            // Loud, per pass: this is the mesh REFUSING to seal a
            // passage, and the user must see it happen, not infer it
            // from a volume delta three stages later.
            std::cerr << "warning: layers: seal guard dropped " << report.sealedStacks.size()
                      << " stacks (" << report.sealedArea << " m2 of wall) on pass " << passesRun
                      << " -- the smoothed landing crossed the midline of a rib thinner than the smoothing"
                      << " diameter there; those faces stop one layer short of the wall, or keep the offset"
                      << " cut where that position is buried too\n";
        }
        if (!report.thinStacks.empty() || report.orphanThinFaces > 0) {
            // Loud, per pass, for the same reason the seal guard is: a
            // collapsed prism is invisible to checkMesh (it reports
            // "Cell volumes OK") and only shows up as a solver Courant
            // collapse, so the refusal must be announced where it
            // happens.
            std::cerr << "warning: layers: collapsed-prism guard terraced " << report.thinStacks.size()
                      << " stacks (" << report.thinArea << " m2 of wall) on pass " << passesRun
                      << " -- the achieved layer height there was below " << minAchievedHeightFrac()
                      << " of the thickness the dict requested; those faces keep the layers they had"
                      << " reached and stay away from the wall\n";
        }
        thinAll.insert(report.thinStacks.begin(), report.thinStacks.end());
        thinAreaAll += report.thinArea;
        orphanThinAll += report.orphanThinFaces;
        sealedAll.insert(report.sealedStacks.begin(), report.sealedStacks.end());
        sealedAreaAll += report.sealedArea;
        orphanSealedAll += report.orphanSealedFaces;
        if (report.badStacks.empty()) break;
        removed.insert(report.badStacks.begin(), report.badStacks.end());
    }
    if (!sealedAll.empty() || orphanSealedAll > 0) {
        std::cout << "layer seal guard: " << sealedAll.size() << " stacks refused their landing, " << sealedAreaAll
                  << " m2 of wall stopped short of the STL, " << orphanSealedAll << " orphan seam faces\n";
    }
    // Collapsed-prism guard, cumulative over every pass, and the
    // shipped-mesh worst achieved/requested first-layer height. Both are
    // printed unconditionally (the second is the number
    // `min_layer_height_frac_above` in Benchmarks/gates.toml binds on --
    // checkMesh demonstrably will not catch a first layer squeezed to a
    // quarter of its requested thickness, so the suite must).
    if (!thinAll.empty() || orphanThinAll > 0) {
        std::cout << "layer collapsed-prism guard: " << thinAll.size() << " stacks terraced, " << thinAreaAll
                  << " m2 of wall kept at the offset cut, " << orphanThinAll << " orphan seam faces\n";
    }
    std::cout << "minLayerHeightFrac = "
              << (result.stats.minLayerHeightFrac < 1e29 ? result.stats.minLayerHeightFrac : -1.0) << "\n";
    if (layerHeightHist()) {
        std::cout << "layer achieved/requested height histogram (bin width 0.05, last bin >= 1.0):\n";
        for (std::size_t j = 0; j < g_heightHistBins.size(); ++j) {
            long tot = 0;
            for (long v : g_heightHistBins[j]) tot += v;
            if (tot == 0) continue;
            long cum = 0;
            for (std::size_t i = 0; i < 21; ++i) {
                cum += g_heightHistBins[j][i];
                std::cout << "  hist L" << j << " " << std::fixed << std::setprecision(2)
                          << (0.05 * static_cast<double>(i)) << " " << g_heightHistBins[j][i] << "  cumFrac "
                          << (tot > 0 ? static_cast<double>(cum) / static_cast<double>(tot) : 0.0) << "\n";
            }
            std::cout << std::defaultfloat << "  hist L" << j << " total " << tot << "\n";
        }
    }
    if (gateDumpOuter) {
        const char* path = std::getenv("NINJA_LAYER_GATE_DUMP");
        std::FILE* f = std::fopen(path, "w");
        if (f) {
            std::fprintf(f, "patchOrdinal,origin,x,y,z,coreLevel\n");
            for (const CondemnedSite& s : allCondemnedSites) {
                std::fprintf(f, "%d,%d,%.9f,%.9f,%.9f,%d\n", s.patchOrdinal, s.origin, s.x, s.y, s.z, s.coreLevel);
            }
            std::fclose(f);
        } else {
            std::cerr << "warning: layer-quality gate: NINJA_LAYER_GATE_DUMP could not open '" << path
                      << "' for writing\n";
        }
        // Sibling terraced-face dump, `<path>.terraced`
        // -- one row per wall face that terraced on the FINAL (converged)
        // pass only (`report` here still holds that last pass's data;
        // earlier passes' terracing is superseded once the gate's
        // removed-set fixpoint is reached and the mesh is re-marched, so
        // dumping only the final pass matches what actually ships, and
        // matches the final `droppedFaces` stat printed below). Same
        // zero-cost-unset / opt-in convention; separate file (not
        // appended to the condemned-sites CSV) so existing downstream
        // tooling that reads the original 6-column format is unaffected.
        const std::string terracedPath = std::string(path) + ".terraced";
        std::FILE* ft = std::fopen(terracedPath.c_str(), "w");
        if (ft) {
            std::fprintf(ft, "patchOrdinal,origin,x,y,z,coreLevel,achievedLayers,specLayers,reason\n");
            for (const TerracedFaceSite& s : report.terracedFaces) {
                std::fprintf(ft, "%d,%d,%.9f,%.9f,%.9f,%d,%d,%d,%s\n", s.patchOrdinal, s.origin, s.x, s.y, s.z,
                             s.coreLevel, s.achievedLayers, s.specLayers, layerDropReasonName(s.reason));
            }
            std::fclose(ft);
        } else {
            std::cerr << "warning: layer-quality gate: NINJA_LAYER_GATE_DUMP terraced sibling '" << terracedPath
                      << "' could not be opened\n";
        }
    }
    if (!removed.empty()) {
        const double frac = report.totalWallFaces0 > 0
                                ? static_cast<double>(removed.size()) / static_cast<double>(report.totalWallFaces0)
                                : 0.0;
        // Format preserved verbatim (%.3f) -- harness-visible line.
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3) << 100.0 * frac;
        std::cout << "layer-quality gate: removed " << removed.size() << " of " << report.totalWallFaces0
                  << " stacks (" << oss.str() << "%, " << initiallyBad << " initially bad, " << passesRun
                  << " passes)\n";
    }
    if (!report.badStacks.empty()) {
        std::cerr << "warning: layer-quality gate: " << report.badStacks.size() << " stacks still bad after "
                  << passesRun << " passes (pass cap reached)\n";
    }
    if (report.orphanBadCells > 0) {
        std::cerr << "warning: layer-quality gate: " << report.orphanBadCells
                  << " bad cells on seam faces (no stack to remove)\n";
    }

    // --- ZERO-FACE PATCH DISCLOSURE (the sealed-flow-passage
    // post-condition) ------------------------------------------------
    //
    // MEASURED (bm_solver_wfp_interfoam, f6e3d05): bm_layers_wfp's
    // `inletB` patch -- the 0.9 x 0.15 m fish-passage inlet slot, the
    // one genuinely thin flow-carrying feature the case dict calls out
    // by name -- ships with ZERO faces and therefore cannot carry its
    // 1.0 m3/s. The same mesh without layers keeps it (17
    // faces). A patch silently going to zero faces is always either a
    // bug or a modelling change the user must know about, and NOTHING
    // in the pipeline said so: the mesh is valid, checkMesh is happy,
    // and the flow domain is simply a different one.
    //
    // The mechanism is NOT the seal guard's (the smoothed landing
    // crossing into the solid). It is arithmetic in the CUT, one stage
    // earlier: cutting the core at d = t thickens every solid by t per
    // side, so a passage narrower than 2t is gone before a prism
    // exists. The slot is 0.15 m; 2t = 0.262 m. That is why the seal
    // guard -- which watches landings, i.e. things that only exist
    // where there is a front to march -- cannot see it, and why the
    // march cannot give it back: by the time this module is called the
    // slot is not in its input mesh either. Hence the disclosure is
    // written against the INPUT as well as the output, so it names the
    // stage responsible instead of leaving the user to guess.
    //
    // Error level, not warning level, and separated from the legitimate
    // case: a DOMAIN-SIDE patch (declared in the dict's `patches` block)
    // is routinely covered by the geometry and empty by design -- this
    // very case declares three of them that way. A GEOMETRY patch is
    // the surface the user named and asked to be meshed; zero faces
    // there is never a design intent.
    try {
        std::set<std::string> domainSideNames;
        for (const PatchSpec& ps : cfg.patches) domainSideNames.insert(ps.name);
        // Counted on the SURVIVING fluid only. `dropDisconnectedCells`
        // runs AFTER this module (main.cpp's pipeline tail), so a patch
        // whose faces all sit on cells the final flood fill is about to
        // discard still has a nonzero nFaces here and would otherwise be
        // missed -- which is exactly how bm_layers_wfp's `inletB` got
        // out: the offset cut narrows the slot's throat until the fluid
        // behind it stops being reachable from locationInMesh, and the
        // patch dies with the component, not with the cut. Same
        // predicate the flood fill itself uses, so this is a prediction
        // of the written mesh rather than an approximation of it.
        const std::vector<char> reachedOut = reachableCellsFromLocation(result.mesh, locationInMesh);
        const std::vector<char> reachedIn = reachableCellsFromLocation(cutMeshIn, locationInMesh);
        auto liveFacesByName = [](const GeneratedMesh& m, const std::vector<char>& reached) {
            std::map<std::string, int> n;
            for (const PatchInfo& p : m.patches) n[p.name] += 0;
            for (int f = m.nInternalFaces; f < m.nFaces(); ++f) {
                const int pid = m.faces.patchId[static_cast<std::size_t>(f)];
                if (pid < 0 || pid >= static_cast<int>(m.patches.size())) continue;
                const int o = m.faces.owner[static_cast<std::size_t>(f)];
                if (o < 0 || !reached[static_cast<std::size_t>(o)]) continue;
                ++n[m.patches[static_cast<std::size_t>(pid)].name];
            }
            return n;
        };
        const std::map<std::string, int> facesIn = liveFacesByName(cutMeshIn, reachedIn);
        const std::map<std::string, int> facesOut = liveFacesByName(result.mesh, reachedOut);
        long zeroGeom = 0, zeroSide = 0;
        for (const auto& kv : facesOut) {
            if (kv.second > 0) continue;
            const std::string& name = kv.first;
            const bool side = domainSideNames.count(name) > 0;
            const auto itIn = facesIn.find(name);
            const int before = itIn != facesIn.end() ? itIn->second : 0;
            std::string type = "?";
            for (const PatchInfo& pi : result.mesh.patches) {
                if (pi.name == name) { type = pi.type; break; }
            }
            if (side) {
                ++zeroSide;
                continue;
            }
            ++zeroGeom;
            std::cerr << "error: layers: patch '" << name << "' (type " << type
                      << ") has ZERO faces on the fluid this run will write"
                      << (before > 0 ? " -- the layer stage lost all " + std::to_string(before)
                                           + " face(s) it was given, by extruding over them or by cutting off"
                                             " the fluid behind them from locationInMesh"
                                     : " -- it already had none on this stage's INPUT, i.e. the offset cut"
                                       " (d = t) closed it: every passage narrower than 2*t is gone before a"
                                       " prism exists")
                      << ". That patch can carry no boundary condition and no flow; if it is meant to, reduce"
                         " the total layer thickness for the surfaces bounding it, or refine locally.\n";
        }
        // Machine-readable for Benchmarks/summary.py (gated per case as
        // zero_face_geometry_patches_below): a patch count, not an area,
        // because losing a patch is a topological fact and a band would
        // hide it.
        std::cout << "zeroFaceGeometryPatches = " << zeroGeom << "\n";
        std::cout << "zeroFaceDomainSidePatches = " << zeroSide << "\n";
    } catch (const std::exception& e) {
        // Disclosure must never be the thing that kills a run that
        // otherwise produced a mesh; say why it could not be computed.
        std::cerr << "warning: layers: zero-face patch disclosure could not be computed: " << e.what() << "\n";
    }
    return result;
}

} // namespace ninja
