// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <vector>

#include "BaseMesh.hpp"
#include "CutData.hpp"
#include "Geometry.hpp"
#include "Stl.hpp"

namespace ninja {

// The layer-by-layer march + prism emission. Consumes the ALREADY-CUT
// mesh (offset wall lands at d = t) and re-homes each layered STL's
// wall patch down to
// the true surface in nLayers graded steps, emitting one prism cell per
// layer-top face per step. The writer never sees
// this module -- it only ever gets a plain GeneratedMesh.
//
// --- GRADING CONTRACT ---
// Geometric division. For nLayers = n and expansionRatio = r > 0, the
// layer thicknesses are
//     layerThickness[j] = t0 * r^j,  j = 0 (AT THE WALL) .. n-1,
//     t0 = t / sum_{j=0}^{n-1} r^j   (so the schedule sums to t exactly),
// i.e. the THINNEST layer is at the wall and each successive layer is r
// times the previous one. r == 1 degenerates to n equal layers; n == 1
// to a single full-thickness step. The march runs OUTERMOST-FIRST
// (step 0 places the outermost interface, the last step lands on the
// wall), and step lengths are expressed as a fraction of each point's
// own REMAINING distance-to-STL rather than as absolute lengths.
//
// Cross-step behaviour (deliberate, and kept):
//   * A face DROPPED at step k is reverted to a plain wall face and
//     rejoins the wall bucket, so step k+1 re-marches it with whatever
//     distance it has left. The result is locally fewer, thicker layers
//     rather than a hole -- graceful degradation, not failure. Such a
//     face is counted ONCE in LayerStats::droppedFaces, however many
//     steps it fails on.
//   * A point FROZEN by the non-inversion clamp at step k re-aims on
//     step k+1 with its remaining d, for the same reason, and likewise
//     counts once.
//   * Only the LAST step snaps to the exact closest point (the landing);
//     intermediate fronts are meant to sit at layer interfaces.
// Measured by the offset_graded gate (Tests/scripts/layer_grading.py).

// A test-only forced-drop sphere (layersDebug{forceDropSphere}):
// any layer-top face whose centroid lies within `radius` of `centre` is
// treated as a validation failure regardless of its geometric validity.
// Test scaffolding, never used outside tests.
struct ForceDropSphere {
    Vec3 centre;
    double radius = 0.0;
};

// One layered STL's march input: `stlIndex` indexes `geom.stls` /
// `Triangle::solidId` / the offset thickness vector already used by
// CutData's offset path; `thickness` is that STL's total t (> 0, or the
// STL would not be in this list at all -- t == 0 STLs never march,
// exactly the L1 "cuts at the true wall" case, still present as an
// ordinary wall patch, untouched by this module).
struct LayerStlSpec {
    int stlIndex = 0;
    std::string wallPatchName;
    double thickness = 0.0;
    // Layer-by-layer march: the march
    // advances in nLayers graded steps (expansionRatio between
    // successive layer thicknesses, thinnest layer at the wall),
    // NEVER the full thickness in one step; each step emits one prism
    // sheet. nLayers 1 degenerates to a single full-thickness step.
    int nLayers = 1;
    double expansionRatio = 1.0;
    // Per-STL smoothing radius fraction (r = smoothRadius * h_local),
    // a first-class dict parameter. Default 1.0: the validated production
    // value -- near-complete marches at minimal sealed area. 0 disables smoothing. Values comparable
    // to a feature's size (in h_local units) can seal narrow passages;
    // keep r well under the smallest feature radius.
    double smoothRadius = 1.0;
};

// Why a wall face ended up without (all of) its layers -- the first refusal it
// met. 1-10 are prismValid's failure codes (Layers.cpp), in its numbering.
enum LayerDropReason : int {
    kDropDegenerateLoop = 1,
    kDropStepTooShort = 2, // a corner marched less than the minimum height
    kDropSideTwist = 3,
    kDropBottomNonConvex = 4,
    kDropVolume = 5,
    kDropBottomFolded = 6,
    kDropSideInverted = 7,
    kDropMeanHeight = 8, // sheared flat
    kDropAspect = 9,
    kDropCollapsed = 10, // achieved height far below the requested step
    kDropSealed = 11,    // seal guard: the landing would bury the face
    kDropForced = 12,    // layersDebug forceDropSphere
    kDropPyramid = 13,   // in-march quality guard, by checkMesh metric
    kDropTet = 14,
    kDropSkew = 15,
    kDropNonOrth = 16,
    kDropCondemned = 17, // its stack was condemned by an earlier gate pass
    kDropNeighbour = 18, // one-ring dilation of a failing face
    kDropHeldSeam = 19,  // shares a spent point with a terrace; never extruded again
    kNumLayerDropReasons = 20
};
inline const char* layerDropReasonName(int r) {
    static const char* const names[kNumLayerDropReasons] = {
        "unknown", "degenerateLoop", "stepTooShort", "sideTwist", "bottomNonConvex", "volume", "bottomFolded",
        "sideInverted", "meanHeight", "aspect", "collapsed", "sealed", "forced", "pyramid",
        "tet", "skew", "nonOrth", "condemned", "neighbour", "heldSeam"};
    return r >= 0 && r < kNumLayerDropReasons ? names[r] : "unknown";
}

// Residual + bookkeeping metrics, reported by --cut-stats and
// normal stdout when `layers` is present.
struct LayerStats {
    double residualMean = 0.0; // area-weighted mean landed distance-to-STL over final wall faces
    double residualMax = 0.0;  // area-weighted... (max is just max, "area-weighted" describes the metric family)
    // DISTINCT counts. A face reverted at step k rejoins the wall
    // bucket and is re-marched at step k+1; counting it again on each
    // failing step would inflate the totals. Both counters key on an
    // ORIGIN id carried in parallel with the wall bucket / front point
    // set, so each originally-distinct wall face (resp. front point) is
    // counted at most once no matter how many steps it fails on.
    int frozenPoints = 0;      // DISTINCT front points the non-inversion clamp bisected down to ~no movement
    int droppedFaces = 0;      // DISTINCT layer-top faces reverted to plain wall faces (drop-locally)
    // Layer-feasibility visibility: summed area of the same DISTINCT
    // `droppedFaces` set, i.e. the wall area where layers were
    // infeasible and the offset cut/seal was kept instead -- the trade
    // is made visible, never silent. Reuses the existing drop
    // accounting -- no separate feasibility mechanism. Both count the
    // fluid reachable from locationInMesh only: wall faces of a component
    // the offset cut disconnected are discarded with it, not dropped.
    double infeasibleWallArea = 0.0;
    // (point, step) pairs the surface guard landed on the nearest triangle
    // because the step crossed it and parity still said fluid beyond (see the
    // guard in Layers.cpp): a handful at thin plates, thousands at a
    // duplicated CAD face.
    long surfaceClampedSteps = 0;
    // The same DISTINCT `droppedFaces` set, split by the FIRST thing that
    // refused the face (index = LayerDropReason). Sums to droppedFaces.
    std::vector<int> droppedByReason;
    // Faces reverted by the in-march quality guard (a prism that would
    // have failed checkMesh's pyramid/tet tests is never emitted).
    // Raw per-pass count, one per (face, step) drop event.
    int qualityDroppedFaces = 0;
    // Per-step observability (raw, NOT de-duplicated): index j is march
    // step j (outermost layer first), summed over all layered STLs.
    std::vector<int> perStepDropped;    // faces reverted at that step
    std::vector<int> perStepFrozen;     // front points frozen at that step
    std::vector<int> perStepPrismCells; // prism cells emitted at that step

    // --- Smoothed-field metrics, all zero when smoothRadius is 0 -- the
    // wired path pays these accumulators only when it is itself active. ---
    //
    // Mechanism-claim numbers: how often the direction
    // consumer had to fall back because the smoothed gradient estimate
    // was degenerate (|grad phi_r| ~ 0) -- falling back to the exact
    // closest-point direction anywhere would silently reintroduce the
    // discontinuity smoothing exists to remove, so this must be
    // counted, never silently substituted.
    long smoothGradFallback = 0;
    // What the smoothing destroys must be loud: count/area of LANDED
    // wall faces where |phi_r - phi| > 0.1 * h_local
    // at at least one of the face's points, plus the area-weighted
    // mean/max of |phi_r - phi| in the SAME (absolute-length) units --
    // reported as a fraction of h_local by the caller, which knows the
    // per-face h_local this struct does not carry.
    long smoothMovedFaces = 0;
    double smoothMovedArea = 0.0;
    double smoothResidWeightedSum = 0.0; // sum over ALL landed wall faces of area * |phi_r - phi|
    double smoothResidAreaSum = 0.0;     // total landed wall area considered (denominator of the above)
    double smoothResidMaxOverH = 0.0;    // max, over all landed points, of |phi_r - phi| / h_local
    // Sealed-region detection: a landed wall face where the smoothed
    // and exact fields DISAGREE ON SIGN (phi_r > 0 but phi < 0, or vice
    // versa) AND the exact field puts that point more than 0.1 * h_local
    // INSIDE the solid -- the smoothing changed inside/outside there by a
    // geometrically meaningful depth, which is where a real passage seals
    // or a rib evaporates -- loud about the geometry, whatever its size
    // in metres. Faces are unioned into connected regions via shared
    // front-point adjacency; one record per region.
    //
    // The depth gate is NOT decoration. A landed wall face sits ON the
    // surface by construction, so phiExact there is ~0 and its sign is
    // numerical noise; an ungated sign test fires on ordinary healthy
    // landings and grossly overstates real sealing, making the metric
    // useless for the radius decision it exists to inform. See
    // Layers.cpp for why |phiExact| is the right gate.
    struct SealedRegion {
        int count = 0;       // faces in the region
        double area = 0.0;   // total area
        Vec3 centroid{0, 0, 0}; // area-weighted world centroid
    };
    std::vector<SealedRegion> smoothSealedRegions;

    // Seal guard: landed wall faces the march REFUSED to
    // emit because their landed points sat inside the solid by more
    // than 0.1 * h_local -- the facing-wall crossing that would have
    // sealed a passage. Same predicate as `smoothSealedRegions` above;
    // these counters record the ACTION, the regions record the
    // measurement. Per-pass (the gate re-marches); the cumulative
    // figure is printed by applyLayers.
    long sealDroppedFaces = 0;
    double sealDroppedArea = 0.0;

    // Collapsed-prism guard (achieved vs REQUESTED layer height).
    // Landed prisms the march refused because their ACHIEVED mean
    // height (volume / top area) fell below kMinAchievedHeightFrac of
    // the layer step the dict asked for -- a cell checkMesh calls
    // "volumes OK" and a VoF solver calls a Courant collapse. Same
    // action as every other guard here: stack condemned, wall keeps the
    // offset cut. Per-pass; applyLayers prints the cumulative figure.
    long thinDroppedFaces = 0;
    double thinDroppedArea = 0.0;
    // The shipped-mesh number the gate binds on: the WORST
    // achieved/requested first-layer height over every wall-layer prism
    // actually emitted. 1.0 means every first layer came out at the
    // thickness the dict asked for; it can exceed 1 where the offset
    // cut left more distance than t. Stays at its sentinel when no
    // wall-layer prism was emitted at all.
    double minLayerHeightFrac = 1e30;

    // --- Landing on phi_r=0, all zero/empty when smoothing is
    // off. ---
    //
    // The phi_r=0 bisection failing to bracket (should not happen
    // -- phi_r > 0 at the pre-landing position and phi_r < 0 inside the
    // solid -- counted rather than silently falling back to the exact
    // snap, so a fallback never reintroduces the exact-landing error
    // invisibly).
    long smoothLandingBisectFallback = 0; // == the two counters below, summed
    long smoothLandingDegenerateGradFallback = 0; // |grad phi_r| ~ 0 at the pre-landing sample
    long smoothLandingNoBracketFallback = 0;       // had a direction but never crossed zero within the cap
    // Landing cost (Illinois-safeguarded regula falsi): total field
    // evaluations spent per candidate (bracket-doubling + regula-falsi
    // iterations, including the initial phi0 sample), summed/maxed/
    // counted across every landing candidate that successfully
    // bracketed -- reported so its typical value (mean = sum/count) is
    // on the record rather than assumed.
    long smoothLandingEvalSum = 0;
    long smoothLandingEvalMax = 0;
    long smoothLandingCandidateCount = 0;
    // Achieved per-stack thickness (origin's step-0 top position,
    // i.e. the exact d=t cut position, to this stack's FINAL landed
    // bottom), one entry per wall-layer front point, collected only at
    // the last march step when smoothing is active. The cut stays
    // at exact d=t while the wall now moves to phi_r=0, so near a
    // convex corner the achieved thickness can fall well short of the
    // nominal t -- a confound worth measuring, not a defect.
    std::vector<double> stackThickness;
    // Threshold counts against EACH point's own spec's nominal t
    // (computed at collection time, since `t` varies per layered STL,
    // unlike `stackThickness` itself which is unit-comparable raw
    // distance and needs no per-spec context downstream).
    long stackThicknessBelowHalfT = 0;
    long stackThicknessBelowQuarterT = 0;
};

// Applies the march + prism emission to `mesh` in place of its layered
// STLs' wall patches. `cfg` supplies the domain bounds (for the
// domain-boundary in-plane projection constraint) and patch names (for
// side faces that land on a domain plane); `perStlBins` is one
// TriangleAabbBins per STL (same ones the offset cut used -- built by
// the caller over `perStlTris`, index-aligned with `geom` STL
// declaration order); `specs` lists only the STLs that actually have a
// `layers` entry (thickness > 0). `cellLevelIn`/`pointLevelIn` are the
// pre-layers level fields (already propagated through the cut); the
// result carries extended level fields for the new prism cells/points.
struct LayersResult {
    GeneratedMesh mesh;
    std::vector<int> cellLevel;
    std::vector<int> pointLevel;
    LayerStats stats;
    // Self-diagnosing wall-distance export: per-point diagnostic
    // flags, aligned 1:1 with
    // `mesh.points` (same size, same indexing) -- so an artifact node in
    // wallDist_<patch>.vtk can be attributed to a specific march-time
    // mechanism instead of guessed at. Empty (size 0) is a valid "not
    // populated" state the caller must check against `mesh.points.size()`
    // before indexing (mirrors every other optional-output convention in
    // this project, e.g. Cutter.hpp's `originCell`).
    // Stack topology export (opt-in): per-CELL stack identity for
    // the final mesh -- `cellStackId` is the originating wall seed face's
    // `faceOrigin` (a "stack" is everything one seed face extrudes over the
    // whole march) and `cellLayerIndex` is the march step that created the
    // cell. Both are -1 for cells that are not layer prisms (core, merged,
    // base). Written out as polyMesh cell fields so a single stack can be
    // isolated in ParaView (threshold on cellStackId) and its ACTUAL topology
    // inspected cell by cell, rather than inferred from aggregate statistics.
    std::vector<int> cellStackId;
    std::vector<int> cellLayerIndex;
    std::vector<char> pointFrozen;  // this point's own origin was ever non-inversion-frozen
    std::vector<char> pointDropped; // this point sits on a face whose origin was ever drop-reverted
};

// `locationInMesh`: the smoothed-field path needs it to sign
// the field via `classifyVertices` (the batched sign source, see
// Geometry.hpp); unused when smoothing is off, so callers pass the
// same locationInMesh they already have (`GeometryConfig::locationInMesh`).
LayersResult applyLayers(const GeneratedMesh& cutMeshIn, const std::vector<int>& cellLevelIn,
                          const std::vector<int>& pointLevelIn, const MeshConfig& cfg,
                          const std::vector<LayerStlSpec>& specs,
                          const std::vector<std::vector<Triangle>>& perStlTris,
                          const std::vector<TriangleAabbBins>& perStlBins,
                          const std::vector<ForceDropSphere>& forceDropSpheres, const Vec3& locationInMesh);

} // namespace ninja
