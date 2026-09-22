// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#include "CutData.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <set>
#include <stdexcept>

#include "Geometry.hpp"

namespace ninja {

namespace {

constexpr double kEps = 1e-9;

// Fixed, ordered fallback ray directions. The primary direction is +x; degenerate hits re-cast with the next
// entry, deterministically, so the choice is a pure function of the
// vertex position + STL (no randomness, no per-cell variation).
const std::vector<Vec3>& fallbackDirections() {
    static const std::vector<Vec3> dirs = [] {
        std::vector<Vec3> raw = {
            Vec3{1.0, 0.0, 0.0},
            Vec3{0.0, 1.0, 0.0},
            Vec3{0.0, 0.0, 1.0},
            Vec3{1.0, 1.0, 1.0},
            Vec3{1.0, -1.0, -1.0},
            Vec3{-1.0, 1.0, -1.0},
            // Deliberately asymmetric directions to avoid systematic
            // grazing hits on symmetric STLs (e.g. an icosphere's
            // triangulation is itself highly axis/diagonal-symmetric,
            // so grid points on symmetry planes can defeat the more
            // "obvious" directions above for every one of them).
            Vec3{0.8471, 0.3941, -0.3585},
            Vec3{-0.2231, 0.9142, 0.3389},
            Vec3{0.4467, -0.7791, 0.4396},
            Vec3{0.6123, 0.1187, 0.7818},
            Vec3{-0.5544, -0.6612, 0.5063},
            Vec3{0.9012, -0.2233, 0.3712},
        };
        for (Vec3& d : raw) {
            const double n = norm(d);
            d = Vec3{d.x / n, d.y / n, d.z / n};
        }
        return raw;
    }();
    return dirs;
}

// Möller–Trumbore ray-triangle intersection. Returns true and sets `t`
// (>= 0 means an intersection along the ray, possibly beyond a finite
// segment; caller checks range) if a valid, non-degenerate intersection
// exists. Sets `degenerate` if the ray grazes an edge/vertex or lies in
// the triangle's plane (barycentric coordinate within kEps of a
// boundary, or the ray near-parallel to the triangle).
bool rayTriangle(const Vec3& origin, const Vec3& dir, const Triangle& tri, double& t, bool& degenerate) {
    degenerate = false;
    const Vec3 e1 = tri.v1 - tri.v0;
    const Vec3 e2 = tri.v2 - tri.v0;
    const Vec3 pvec = cross(dir, e2);
    const double det = dot(e1, pvec);
    // Relative parallel test -- det scales with |dir||e1||e2| (see the
    // matching comment in computeEdgeIntercepts; `dir` is unit here but
    // real-CAD triangle edges can be tiny, so the product still matters).
    if (std::fabs(det) < kEps * norm(dir) * norm(e1) * norm(e2)) {
        return false; // parallel; not degenerate unless ray lies in-plane and would hit -
                       // treated conservatively below by the caller's u/v/w check on other tris.
    }
    const double invDet = 1.0 / det;
    const Vec3 tvec = origin - tri.v0;
    const double u = dot(tvec, pvec) * invDet;
    const Vec3 qvec = cross(tvec, e1);
    const double v = dot(dir, qvec) * invDet;
    const double w = 1.0 - u - v;
    if (u < -kEps || v < -kEps || w < -kEps || u > 1.0 + kEps || v > 1.0 + kEps || w > 1.0 + kEps) {
        return false; // clean miss
    }
    t = dot(e2, qvec) * invDet;
    // Near a triangle edge or vertex: degenerate, caller must re-cast —
    // UNLESS the near-hit is at (or behind) the ray origin itself
    // (t <= 0): a classification point that coincidentally sits exactly
    // on an STL vertex/edge always looks "degenerate" at t=0 for every
    // incident triangle regardless of direction (every fallback would
    // exhaust); such a t<=0 graze contributes no crossing either way,
    // so it is safe to treat as a clean non-crossing rather than force
    // a retry.
    if ((std::fabs(u) < kEps || std::fabs(v) < kEps || std::fabs(w) < kEps) && t > kEps) {
        degenerate = true;
        return false;
    }
    return t > kEps;
}

// Ray-parity crossing count from `origin` in direction `dir`, over the
// candidate triangle indices `idxs` (a superset of every triangle the
// ray could geometrically hit, per queryRayTriangles' binning -- see
// Geometry.hpp). Returns {count, ok}; ok=false means a degenerate hit
// occurred and the caller should retry with the next direction.
std::pair<int, bool> parityCount(const Vec3& origin, const Vec3& dir, const std::vector<Triangle>& tris,
                                  const std::vector<int>& idxs) {
    int count = 0;
    for (int i : idxs) {
        const Triangle& tri = tris[static_cast<std::size_t>(i)];
        double t = 0.0;
        bool degenerate = false;
        const bool hit = rayTriangle(origin, dir, tri, t, degenerate);
        if (degenerate) {
            return {0, false};
        }
        if (hit && t > kEps) {
            ++count;
        }
    }
    return {count, true};
}

// `bins` is an AABB-based accelerator over the SAME `tris` (see
// Geometry.hpp's TriangleAabbBins/queryRayTriangles): for each fixed
// fallback direction, only the triangles whose bins the ray's path
// actually crosses are tested, instead of the full soup. Exactness:
// queryRayTriangles returns a superset of every triangle the ray could
// hit (see its docstring), so `count` here is bit-identical to the old
// full-soup scan -- WHICH triangles get tested changes, HOW they get
// tested (rayTriangle, epsilons, tie-breaking) does not.
//
// Multi-ray parity voting -- LONE-OUTLIER override of the first ray's
// verdict, nothing more.
//
// A ray that is NOT flagged degenerate (no edge/vertex graze --
// rayTriangle's own `degenerate` stays false) can still return a
// spuriously WRONG parity on open (non-watertight) CAD when it
// travels exactly along a genuine surface feature (e.g. a
// tire/ground contact line a few 1e-5 above the ray): it passes
// straight through the open shells and racks up an ODD count of
// real, non-degenerate crossings, flooding a point far from any
// triangle as "solid". A "first non-degenerate direction wins" rule
// cannot detect this -- the ray is not degenerate by this function's
// (correct, narrower) definition, it is just wrong because the input
// geometry is open. Rejecting-and-recasting near-tolerance hits
// would not help either: such a ray does not graze an edge, it runs
// through actual holes in the surface, arbitrarily far from any
// triangle boundary.
//
// So: cast EVERY fixed fallback direction (still a pure,
// deterministic function of point + STL -- no randomness, no
// per-cell/per-rank variation, so this stays np-invariant by
// construction) and keep the FIRST non-degenerate ray's verdict
// UNLESS that verdict is a strict LONE OUTLIER: no other
// non-degenerate ray agrees with it, and at least two other rays
// unanimously contradict it. Only then flip to the unanimous
// opposition.
//
// Why not a plain majority (or a closest-point-normal tie-break)?
// Near a dirty-CAD surface the votes can genuinely split (7 vs 5 was
// observed): ray parity is structurally ambiguous there and NO
// aggregate of it is more trustworthy than any other -- and flipping
// such points away from the first-ray answer reintroduces the very
// "cut edge reports no STL crossing" fatal this rule exists to
// remove. So any point whose first verdict has even ONE supporting
// ray keeps it.
//
// The rule adds no epsilon knob (the trigger is exact integer vote
// counts: firstVerdict's side == 1 vote, opposition >= 2 votes,
// unanimously) and is a no-op by construction on any STL where all
// rays agree -- i.e. every watertight input.
bool isSolidAt(const Vec3& p, const std::vector<Triangle>& tris, const TriangleAabbBins& bins) {
    int solidVotes = 0;
    int fluidVotes = 0;
    bool haveFirst = false;
    bool firstSolid = false;
    for (const Vec3& d : fallbackDirections()) {
        const std::vector<int> candidates = queryRayTriangles(bins, p, d);
        auto [count, ok] = parityCount(p, d, tris, candidates);
        if (!ok) {
            continue; // degenerate ray: no vote, try the next direction
        }
        const bool thisSolid = (count % 2) == 1;
        if (!haveFirst) {
            haveFirst = true;
            firstSolid = thisSolid;
        }
        if (thisSolid) {
            ++solidVotes;
        } else {
            ++fluidVotes;
        }
        // Early exit (pure optimization, provably outcome-identical):
        // the lone-outlier flip below requires the first verdict to
        // have NO other supporter -- as soon as a second ray agrees
        // with the first, the answer is fixed at `firstSolid`
        // regardless of every remaining ray, so stop casting. Typical
        // (clean-geometry) cost is therefore 2 rays, not the full list;
        // the old code's cost was 1.
        if ((firstSolid ? solidVotes : fluidVotes) >= 2) {
            return firstSolid;
        }
    }
    if (!haveFirst) {
        throw std::runtime_error(
            "CutData error: all fallback ray directions were degenerate for classification point (" +
            std::to_string(p.x) + " " + std::to_string(p.y) + " " + std::to_string(p.z) + ")");
    }
    // Lone-outlier override (see the function comment above for the
    // full rationale): flip the first ray's verdict only when NO other
    // non-degenerate ray agrees with it and at least TWO rays
    // unanimously contradict it. Every other configuration (unanimous
    // agreement, any support at all for the first verdict, or a lone
    // opposing ray) keeps the first-ray answer exactly.
    const int firstSideVotes = firstSolid ? solidVotes : fluidVotes;
    const int otherSideVotes = firstSolid ? fluidVotes : solidVotes;
    if (firstSideVotes == 1 && otherSideVotes >= 2) {
        return !firstSolid;
    }
    return firstSolid;
}

// True iff `p` lies within kEps of the plane of `tri` AND within its
// triangle footprint (barycentric coords all >= -kEps). Same
// plane+footprint predicate already used by computeEdgeIntercepts's
// degenerate-edge fallback; hoisted here so both the
// intercept-snap path and vertex classification share one definition of
// "on the STL surface" (single-sourced, deterministic).
bool pointOnTriangle(const Vec3& p, const Triangle& tri) {
    const Vec3 e1 = tri.v1 - tri.v0;
    const Vec3 e2 = tri.v2 - tri.v0;
    const Vec3 n = cross(e1, e2);
    const double nlen2 = dot(n, n);
    if (nlen2 < kEps) {
        return false; // degenerate triangle
    }
    const Vec3 w = p - tri.v0;
    if (std::fabs(dot(w, n)) > kEps * std::sqrt(nlen2)) {
        return false; // not on the triangle's plane
    }
    const double d00 = dot(e1, e1);
    const double d01 = dot(e1, e2);
    const double d11 = dot(e2, e2);
    const double d20 = dot(w, e1);
    const double d21 = dot(w, e2);
    const double denom = d00 * d11 - d01 * d01;
    if (std::fabs(denom) < kEps) {
        return false;
    }
    const double v = (d11 * d20 - d01 * d21) / denom;
    const double ww = (d00 * d21 - d01 * d20) / denom;
    const double u = 1.0 - v - ww;
    return !(u < -kEps || v < -kEps || ww < -kEps);
}

// `bins` restricts the scan to triangles whose (padded) AABB contains
// (a tiny box around) `p` -- a superset of every triangle
// pointOnTriangle could possibly accept: pointOnTriangle's own
// tolerances (plane distance kEps*|n|, barycentric slop kEps on
// dimensionless u/v/w) are bounded by kEps (1e-9) times the triangle's
// own edge lengths/area, which is many orders of magnitude smaller than
// the bin grid's numerical padding (~cellSize*1e-6) for any STL this
// project's bin sizing produces -- see Geometry.hpp's
// queryAabbTriangles docstring for the general AABB-overlap superset
// argument.
// Input-precision on-surface band: binary STLs store float32
// coordinates, so a wall plane at an exact grid-vertex coordinate can be
// stored as the nearest float32, up to ~1e-7 * |coordinate| away from
// the true (double) grid vertex -- far above kEps (1e-9). A vertex
// within this band of the surface cannot be meaningfully classified by
// parity (the STL doesn't resolve position more finely than its own
// storage precision), so it is treated as solid by the same "on-surface
// is solid" convention as the exact case, just widened to the input's
// real precision (see classifyVertices below for the full rationale).
// Single-sourced here so classifyVertices' vertex classification and
// computeEdgeIntercepts' degenerate-edge snap fallback agree on exactly
// the same band -- a vertex classified solid BY the band must also be
// recognized as "on the surface" BY the snap fallback, or a cut edge
// whose true STL intersection point falls just outside the segment
// (rounded the "wrong" way by the same float32 noise) reports a false
// "no STL crossing" error instead of snapping to that solid endpoint.
constexpr double kInputPrecisionEps = 4.0 * 1.1920929e-7;
double inputPrecisionBand(const TriangleAabbBins& bins) {
    double maxAbs = 0.0;
    maxAbs = std::max({maxAbs, std::fabs(bins.lo.x), std::fabs(bins.lo.y), std::fabs(bins.lo.z),
                        std::fabs(bins.hi.x), std::fabs(bins.hi.y), std::fabs(bins.hi.z)});
    return kInputPrecisionEps * maxAbs;
}

bool isOnSurface(const Vec3& p, const std::vector<Triangle>& tris, const TriangleAabbBins& bins) {
    for (int i : queryAabbTriangles(bins, p, p)) {
        if (pointOnTriangle(p, tris[static_cast<std::size_t>(i)])) {
            return true;
        }
    }
    return false;
}

} // namespace

std::vector<bool> classifyVertices(const std::vector<Vec3>& points, const std::vector<Triangle>& tris,
                                    const Vec3& locationInMesh) {
    // AABB-based triangle binning accelerator (see Geometry.hpp) --
    // replaces an O(nVertices x nTriangles) brute-force scan with
    // per-query candidate lists that are a provable superset of the
    // full-soup scan, so results are bit-identical. Built once
    // per call (shared by every vertex/ray below), same "single query
    // structure per shared computation" discipline as the rest of
    // CutData.cpp.
    const TriangleAabbBins bins = buildTriangleAabbBins(tris);
    const bool locationSolid = isSolidAt(locationInMesh, tris, bins);
    // See inputPrecisionBand()'s docstring above: a grid vertex within
    // this band of the surface used to classify FLUID (float32 STL
    // noise exceeds the old 1e-9 on-surface tolerance) and the resulting
    // cut edge intercept landed ~band-scale from the kept corner --
    // "Edges too small" in checkMesh. `band` is ~1e-5 * the model's
    // coordinate magnitude here (4 float32 machine-epsilons), negligible
    // next to the smoothing radius (~h) that also consumes this same
    // classification as its sign source.
    const double band = inputPrecisionBand(bins);
    // Each vertex verdict is a pure function of (position, STL) -- the loop
    // parallelizes with no shared mutable state. `solidBuf` is a char buffer
    // (not vector<bool>) so concurrent element writes never touch shared
    // bits; verdicts are bit-identical to the serial sweep for any schedule.
    std::vector<char> solidBuf(points.size(), 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256)
#endif
    for (std::ptrdiff_t pi = 0; pi < static_cast<std::ptrdiff_t>(points.size()); ++pi) {
        const std::size_t i = static_cast<std::size_t>(pi);
        // Deterministic on-surface convention: a grid vertex lying
        // exactly on (within kEps of) the STL surface has a parity
        // result that depends on which side of the surface the
        // fallback ray direction happens to graze, which is NOT a
        // stable classification -- it can put a whole coincident face
        // of such vertices on the fluid side, silently dropping cut
        // geometry that should be a wall (and net kept-volume error
        // can stay near-zero because the misclassification is locally
        // compensating, not globally absent). Detecting "on surface" BEFORE the parity test and
        // classifying such vertices as SOLID by convention side-steps
        // the ambiguous ray test entirely: a solid vertex on the
        // surface produces cut edges that snap exactly to it (the
        // existing snap-to-endpoint logic in computeEdgeIntercepts), so
        // the wall lands exactly on the coincident surface instead of
        // vanishing. Single-sourced per vertex (evaluated once, shared
        // by every incident cell) and deterministic (pure function of
        // vertex position + STL, no per-cell variation).
        // anyTriangleWithin(..., band) subsumes the old exact
        // isOnSurface check: band is computed from the input's float32
        // noise floor and is never smaller than pointOnTriangle's own
        // kEps=1e-9 tolerance for any STL this project handles, so a
        // single distance check replaces the former two-tier test.
        if (anyTriangleWithin(bins, points[i], band)) {
            solidBuf[i] = 1;
            continue;
        }
        const bool pointSolid = isSolidAt(points[i], tris, bins);
        // Vertex is solid iff it differs from the known-fluid reference
        // point's inside/outside state.
        solidBuf[i] = (pointSolid != locationSolid) ? 1 : 0;
    }
    return std::vector<bool>(solidBuf.begin(), solidBuf.end());
}

std::unordered_map<EdgeKey, Intercept, EdgeKeyHash> computeEdgeIntercepts(
    const std::vector<Vec3>& points, const std::vector<Triangle>& tris,
    const std::vector<bool>& vertexSolid, const std::vector<EdgeKey>& edges) {
    std::unordered_map<EdgeKey, Intercept, EdgeKeyHash> result;
    // AABB-based accelerator (see Geometry.hpp / classifyVertices above);
    // built once and reused for every edge's segment query and the
    // degenerate-tie-break point queries below.
    const TriangleAabbBins bins = buildTriangleAabbBins(tris);
    // Input-precision band, shared by the degenerate fallback's
    // on-surface test and the solid-endpoint intercept snap below --
    // must match classifyVertices' band exactly (single source).
    const double band = inputPrecisionBand(bins);

    for (const EdgeKey& e : edges) {
        const int a = e.first;
        const int b = e.second;
        if (vertexSolid[static_cast<std::size_t>(a)] == vertexSolid[static_cast<std::size_t>(b)]) {
            continue; // uncut edge, even if the STL grazes it (accepted MVP limit)
        }
        const Vec3& pa = points[static_cast<std::size_t>(a)];
        const Vec3& pb = points[static_cast<std::size_t>(b)];
        const Vec3 dir = pb - pa;

        // Endpoint nearest the solid side: prefer crossings close to it
        // (deterministic; keeps the most fluid).
        const bool aSolid = vertexSolid[static_cast<std::size_t>(a)];

        bool found = false;
        double bestT = 0.0;
        int bestSolidId = 0;
        // Candidate triangles whose (padded) AABB overlaps the edge
        // segment's AABB -- a superset of every triangle the brute-force
        // loop below would have tested (see queryAabbTriangles'
        // docstring); the math per candidate is byte-for-byte unchanged.
        const Vec3 segLo{std::min(pa.x, pb.x), std::min(pa.y, pb.y), std::min(pa.z, pb.z)};
        const Vec3 segHi{std::max(pa.x, pb.x), std::max(pa.y, pb.y), std::max(pa.z, pb.z)};
        for (int triIdx : queryAabbTriangles(bins, segLo, segHi)) {
            const Triangle& tri = tris[static_cast<std::size_t>(triIdx)];
            const Vec3 e1 = tri.v1 - tri.v0;
            const Vec3 e2 = tri.v2 - tri.v0;
            const Vec3 pvec = cross(dir, e2);
            const double det = dot(e1, pvec);
            // RELATIVE parallel test: det has units of L^3 (edge x triangle
            // area), so an ABSOLUTE 1e-9 threshold silently rejects
            // every triangle once refined sub-edges (<1mm) meet real
            // CAD-sized triangles -- a genuinely cut edge then reports
            // "no STL crossing". Scale by |dir||e1||e2|.
            const double detScale = norm(dir) * norm(e1) * norm(e2);
            if (std::fabs(det) < kEps * detScale) {
                continue;
            }
            const double invDet = 1.0 / det;
            const Vec3 tvec = pa - tri.v0;
            const double u = dot(tvec, pvec) * invDet;
            if (u < -kEps || u > 1.0 + kEps) {
                continue;
            }
            const Vec3 qvec = cross(tvec, e1);
            const double v = dot(dir, qvec) * invDet;
            const double w = 1.0 - u - v;
            if (v < -kEps || w < -kEps || v > 1.0 + kEps || w > 1.0 + kEps) {
                continue;
            }
            const double t = dot(e2, qvec) * invDet;
            if (t < -kEps || t > 1.0 + kEps) {
                continue;
            }
            const double tc = std::clamp(t, 0.0, 1.0);
            // Nearest-to-solid-endpoint: minimize t if a is solid, else
            // maximize t (closest to b, the solid endpoint).
            if (!found) {
                found = true;
                bestT = tc;
                bestSolidId = tri.solidId;
            } else if (aSolid ? (tc < bestT) : (tc > bestT)) {
                bestT = tc;
                bestSolidId = tri.solidId;
            }
        }

        if (!found) {
            // Degenerate tie-break fallback: when an STL vertex
            // or a whole STL face sits bit-exact on/along a grid
            // vertex/edge, the entire edge (or the coplanar portion of
            // it) is numerically degenerate for the standard
            // Moller-Trumbore test above (det==0 -- the ray direction
            // lies IN the triangle's plane), even though the edge is
            // genuinely cut (endpoint parities differ) and one endpoint
            // is itself a point ON the STL surface -- OR within
            // classifyVertices' same input-precision band of it (the
            // vertex that classified solid via the band, not via an
            // exact match: its true nearest triangle can lie just
            // outside the edge segment, e.g. a float32-rounded plate at
            // x=0.300000012 with the segment ending exactly at the
            // double grid vertex x=0.3 -- the Moller-Trumbore scan above
            // then correctly finds no in-segment crossing). This is the
            // same class of residual, per-entity degeneracy as an
            // intercept landing exactly on a grid vertex, which snaps to
            // that vertex; handling it by
            // extending the snap-to-endpoint logic keeps the fix local
            // to one edge. Test: is pa (or pb) within `band` of, and
            // within the (band-padded) footprint of, any triangle?
            // Union of pa's and pb's `band`-expanded candidate
            // triangles, kept in ascending triangle-index order
            // (queryAabbTriangles already returns sorted-unique lists;
            // merge preserves that) so the scan order below exactly
            // matches the old triangle-index-major "for tri: check pa
            // then pb" loop -- required for bit-identical tie-breaking,
            // not just the same final result.
            const Vec3 bandVec{band, band, band};
            std::vector<int> candA = queryAabbTriangles(bins, pa - bandVec, pa + bandVec);
            std::vector<int> candB = queryAabbTriangles(bins, pb - bandVec, pb + bandVec);
            std::vector<int> cand;
            std::set_union(candA.begin(), candA.end(), candB.begin(), candB.end(), std::back_inserter(cand));
            for (int triIdx : cand) {
                const Triangle& tri = tris[static_cast<std::size_t>(triIdx)];
                for (int which = 0; which < 2; ++which) {
                    const Vec3& p = which == 0 ? pa : pb;
                    if (!pointOnTriangle(p, tri) && pointTriDistSq(p, tri) > band * band) {
                        continue;
                    }
                    found = true;
                    bestT = which == 0 ? 0.0 : 1.0;
                    bestSolidId = tri.solidId;
                    break;
                }
                if (found) {
                    break;
                }
            }
        }

        if (!found) {
            // Watertight-STL assumption guarantees a crossing when
            // endpoint statuses differ; this should not happen for
            // well-behaved input (no STL validation is performed).
            throw std::runtime_error(
                "CutData error: cut edge reports no STL crossing (malformed STL?) edge (" +
                std::to_string(pa.x) + " " + std::to_string(pa.y) + " " + std::to_string(pa.z) +
                ") -> (" + std::to_string(pb.x) + " " + std::to_string(pb.y) + " " +
                std::to_string(pb.z) + ")");
        }

        // Snap-to-endpoint: a crossing within machine tolerance of an
        // endpoint keeps the edge "cut" but takes that endpoint's exact
        // coordinates.
        Vec3 pt = pa + dir * bestT;
        if (bestT < kEps) {
            pt = pa;
        } else if (bestT > 1.0 - kEps) {
            pt = pb;
        }
        // Input-precision snap, SOLID endpoint only: when the crossing
        // lies within the float32 noise band of the solid endpoint
        // (a vertex the band classification put ON the surface), every
        // edge incident to that vertex must emit the IDENTICAL intercept
        // position, or pairs of intercepts a float32-ulp apart become
        // checkMesh "Edges too small". Snapping to a solid endpoint is
        // safe -- a solid corner is never a kept mesh point, so no
        // kept-point/intercept coincidence can arise (the failure mode
        // that forbids snapping toward the FLUID endpoint).
        const Vec3& solidP = aSolid ? pa : pb;
        const Vec3 dSolid = pt - solidP;
        if (dot(dSolid, dSolid) < band * band) {
            pt = solidP;
        }
        result[e] = Intercept{pt, bestSolidId};
    }
    return result;
}

namespace {

// Shared per-point predicate for the implicit offset solid
// { p : insideOriginal(p) OR (exists s: distToStl_s(p) <= stlThickness[s]) }
// -- the single source of truth for "is this point in the offset solid",
// evaluated by the half-grid classifier at two radii.
bool offsetSolidAt(const Vec3& p, const std::vector<Triangle>& tris, const TriangleAabbBins& combinedBins,
                    bool locationSolid, const std::vector<TriangleAabbBins>& perStlBins,
                    const std::vector<double>& stlThickness) {
    // On-surface convention (same rationale as classifyVertices' plain
    // path): a point exactly on the true STL surface is SOLID,
    // side-stepping the ambiguous ray test. anyTriangleWithin's `<=`
    // comparison below independently covers the |d - t| coincidence
    // case for t > 0.
    if (isOnSurface(p, tris, combinedBins)) {
        return true;
    }
    const bool pointSolid = isSolidAt(p, tris, combinedBins);
    if (pointSolid != locationSolid) {
        return true;
    }
    for (std::size_t s = 0; s < stlThickness.size(); ++s) {
        if (stlThickness[s] > 0.0 && s < perStlBins.size() && anyTriangleWithin(perStlBins[s], p, stlThickness[s])) {
            return true;
        }
    }
    return false;
}

std::vector<TriangleAabbBins> buildPerStlBins(const std::vector<std::vector<Triangle>>& perStlTris,
                                               const std::vector<double>& stlThickness) {
    std::vector<TriangleAabbBins> perStlBins(perStlTris.size());
    for (std::size_t s = 0; s < perStlTris.size(); ++s) {
        if (s < stlThickness.size() && stlThickness[s] > 0.0) {
            perStlBins[s] = buildTriangleAabbBins(perStlTris[s]);
        }
    }
    return perStlBins;
}

} // namespace

std::vector<bool> classifyVerticesOffsetHalfGrid(const std::vector<Vec3>& points,
                                                 const std::vector<Triangle>& tris,
                                                 const std::vector<std::vector<Triangle>>& perStlTris,
                                                 const Vec3& locationInMesh,
                                                 const std::vector<double>& stlThickness,
                                                 const std::vector<double>& bandPerPoint,
                                                 std::vector<char>& onOut) {
    const TriangleAabbBins combinedBins = buildTriangleAabbBins(tris);
    const std::vector<TriangleAabbBins> perStlBins = buildPerStlBins(perStlTris, stlThickness);
    const bool locationSolid = isSolidAt(locationInMesh, tris, combinedBins);

    // The three states are read off the SAME exact predicate the binary path
    // uses, evaluated at two radii -- no new distance query, no epsilon:
    //   near = offsetSolid at (t - band)   ->  IN
    //   far  = offsetSolid at (t + band)   ->  not OUT
    //   far && !near                       ->  ON  (|phi| < band)
    // Shrunk thicknesses are floored at 0: an STL thinner than the band
    // contributes only its own parity term, which is what t == 0 already means
    // to offsetSolidAt. The band is PER VERTEX (h/4 of the finest leaf
    // touching it, see main.cpp); the bins do not depend on the radius (only
    // on thickness > 0), so the one unshrunk set serves both evaluations.

    // Pure per-vertex verdicts; parallel over points (char buffer, not
    // vector<bool>, so concurrent writes never share bits -- results are
    // bit-identical to the serial sweep for any schedule).
    std::vector<char> solidBuf(points.size(), 0);
    onOut.assign(points.size(), 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256)
#endif
    for (std::ptrdiff_t pi = 0; pi < static_cast<std::ptrdiff_t>(points.size()); ++pi) {
        const std::size_t i = static_cast<std::size_t>(pi);
        const double b = bandPerPoint[i];
        std::vector<double> shrunk(stlThickness.size()), grown(stlThickness.size());
        for (std::size_t s = 0; s < stlThickness.size(); ++s) {
            shrunk[s] = stlThickness[s] > 0.0 ? std::max(0.0, stlThickness[s] - b) : 0.0;
            grown[s] = stlThickness[s] > 0.0 ? stlThickness[s] + b : 0.0;
        }
        const bool nearSolid = offsetSolidAt(points[i], tris, combinedBins, locationSolid, perStlBins, shrunk);
        const bool farSolid = offsetSolidAt(points[i], tris, combinedBins, locationSolid, perStlBins, grown);
        // ON counts as FLUID (kept), not solid. Marking ON as solid is wrong
        // in two compounding ways: it dilates the removed region by the ON
        // band and deletes whole cells at corners with no cut face at all;
        // and if only cell STATUS ignores ON while clipFace still reads the
        // same mask, faces around ON vertices get clipped while the cell
        // calls itself uncut and emits ORIGINAL loops -- one cell emitting a
        // clipped face while its neighbour emits the full one, breaking edge
        // closure.
        // ON-as-fluid is only safe because computeEdgeInterceptsOffsetHalfGrid
        // places an ON intercept AT the vertex and Cutter.cpp remaps it onto
        // that vertex's own point INDEX: the kept corner and the intercept are
        // then literally the same point, and dedupLoop collapses the repeat.
        solidBuf[i] = nearSolid ? 1 : 0;           // IN only
        onOut[i] = (farSolid && !nearSolid) ? 1 : 0;
    }
    return std::vector<bool>(solidBuf.begin(), solidBuf.end());
}

// Half-grid intercepts: no sampling, no bisection. Exactly one endpoint of a
// status-changing edge is solid (ON or IN, both marked in `vertexSolid`).
//   solid endpoint is ON -> the surface passes THROUGH it; the intercept is
//                           that vertex's own position, and Cutter.cpp maps
//                           such an intercept back onto the grid vertex's
//                           existing point index so several edges meeting at
//                           one ON vertex cannot spawn coincident points.
//   solid endpoint is IN -> the crossing is interior to the edge; snap to the
//                           MIDPOINT (the nearest half-grid node).
std::unordered_map<EdgeKey, Intercept, EdgeKeyHash> computeEdgeInterceptsOffsetHalfGrid(
    const std::vector<Vec3>& points, const std::vector<Triangle>& tris, const std::vector<bool>& vertexSolid,
    const std::vector<char>& onSurface, const std::vector<EdgeKey>& edges, OffsetCutStats& stats) {
    std::unordered_map<EdgeKey, Intercept, EdgeKeyHash> result;
    const TriangleAabbBins combinedBins = buildTriangleAabbBins(tris);
    stats.multiRootEdges = 0;
    for (const EdgeKey& e : edges) {
        const std::size_t a = static_cast<std::size_t>(e.first);
        const std::size_t b = static_cast<std::size_t>(e.second);
        if (vertexSolid[a] == vertexSolid[b]) continue;
        // Exactly one endpoint is IN (solid); the other is ON or OUT. If it is
        // ON the surface passes through it -- the intercept IS that vertex.
        // Otherwise the crossing is interior to the edge: snap to the midpoint,
        // the nearest half-grid node.
        const std::size_t fluidEnd = vertexSolid[a] ? b : a;
        const Vec3 pt = onSurface[fluidEnd] ? points[fluidEnd] : (points[a] + points[b]) * 0.5;
        int solidId = 0;
        if (!tris.empty()) {
            const ClosestHit hit = closestPointOnSoup(combinedBins, pt);
            if (hit.triangle >= 0) solidId = tris[static_cast<std::size_t>(hit.triangle)].solidId;
        }
        result[e] = Intercept{pt, solidId};
        ++stats.quantizedIntercepts;
        if (onSurface[fluidEnd]) ++stats.onVertexIntercepts;
    }
    return result;
}


std::vector<EdgeKey> collectEdges(const GeneratedMesh& mesh) {
    std::set<EdgeKey> unique;
    for (int f = 0; f < mesh.nFaces(); ++f) {
        const IntSpan pts = mesh.faces.pointsOf(f);
        const int n = pts.size();
        for (int i = 0; i < n; ++i) {
            const int a = pts[i];
            const int b = pts[(i + 1) % n];
            unique.insert(makeEdgeKey(a, b));
        }
    }
    return std::vector<EdgeKey>(unique.begin(), unique.end());
}

} // namespace ninja
