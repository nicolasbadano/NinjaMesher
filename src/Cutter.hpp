// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <vector>

#include "BaseMesh.hpp"
#include "CutData.hpp"
#include "Stl.hpp"

namespace ninja {

// One smoothing sealed region -- a connected patch of landed wall
// faces where the smoothed and exact fields DISAGREE ON SIGN, i.e. where the
// smoothing moved the surface across itself and a passage sealed or a rib
// evaporated. Carried per region (not just summed) so the radius sweep can ask
// whether the SAME places seal as r grows.
struct SealedRegionRecord {
    int count = 0;
    double area = 0.0;
    Vec3 centroid{0, 0, 0};
};

// --cut-stats numbers: kept/removed/cut cell counts and total kept
// fluid volume, alongside the vertex/edge counts.
struct CutStats {
    int solidVertices = 0;
    int fluidVertices = 0;
    int cutEdges = 0;
    int keptCells = 0;
    int cutCells = 0;
    int removedCells = 0;
    double totalVolume = 0.0;
    double wallArea = 0.0; // summed area of all STL wall-patch faces
    // Cells that failed cutMesh's closed-polyhedron post-condition (< 4
    // faces, an edge of incidence != 2, or non-positive volume). Removed
    // rather than emitted unless NINJA_CORE_VALIDITY=count.
    int invalidCoreCells = 0;
    // Subset of invalidCoreCells: cells whose emitted face set was a closed
    // polyhedron of positive volume but failed the checkMesh-mimicking
    // per-face geometry tests (pyramid orientation about the cell centre,
    // vertex-anchored tet decomposition, relative-area degeneracy). Removed
    // for the same reason -- see cellGeometryIsValid in Cutter.cpp.
    int invalidGeometryCells = 0;
    // Cut cells removed because a face they would ship as a BOUNDARY face
    // breaches checkMesh's boundary skewness -- see the post-pass in cutMesh.
    int skewBoundaryCells = 0;
    // Half-grid cut counters -- see CutData.hpp's OffsetCutStats. Populated
    // by the caller when a half-grid CutData was used; cutMesh itself does
    // not know which path its CutData came from and never touches them.
    int multiRootEdges = 0;
    int quantizedIntercepts = 0;
    int onVertexIntercepts = 0;
    // Count of coplanar flaps re-parented to the cell on the far
    // side of their host face -- see repairCoplanarFlaps. Zero on every mesh
    // where the chord chain never grazes a foreign ON sheet.
    int repairedCoplanarFlaps = 0;
    // Diagnostic tripwire for the invariant the flap pass restores: cut faces
    // whose pyramid volume about their own cell's centroid is <= 0 (checkMesh's
    // "incorrectly oriented face pyramid"). Reported, never repaired -- a
    // non-zero value means a NEW degeneracy family reached the cut, and the
    // repair belongs where the loop is built, not in a post-hoc face flip.
    int badCutFacePyramids = 0;
    // The layer march's residual + bookkeeping metrics
    // (points that never reached their target step get clamped short, and that
    // shortfall/count is what these fields track). Zero unless a `layers`
    // block is present -- populated by main.cpp's runPipeline from
    // Layers.cpp's LayerStats after applyLayers runs. cutMesh itself
    // never touches these (same "populated by the caller" pattern as
    // multiRootEdges above).
    double residualMean = 0.0;
    double residualMax = 0.0;
    int frozenPoints = 0;
    int droppedFaces = 0;
    int qualityDroppedFaces = 0; // in-march quality-guard drops (see LayerStats)
    std::vector<int> droppedByReason; // droppedFaces by first refusal (see LayerStats, LayerDropReason)
    long surfaceClampedSteps = 0;     // see LayerStats
    long buriedLandingFaces = 0;      // see LayerStats
    double buriedLandingArea = 0.0;
    double buriedLandingMaxDepth = 0.0;
    double buriedLandingMaxDepthOverH = 0.0;
    // Per-march-step observability (raw counts, not
    // de-duplicated); frozenPoints/droppedFaces above are DISTINCT
    // counts. Empty unless a `layers` block is present.
    std::vector<int> perStepDropped;
    std::vector<int> perStepFrozen;
    std::vector<int> perStepPrismCells;
    // Sliver merge (offset-cut path only). A sub-threshold
    // cut cell is merged into a kept neighbour instead of removed (see
    // mergeSlivers below); mergeFallbackRemovals counts sliver cells
    // whose connected sliver component had no kept neighbour at all
    // (degenerate case, falls back to the old removal treatment).
    int mergedSliverCells = 0;
    int mergeFallbackRemovals = 0;
    // Geometry-invalid slivers mergeSlivers could not merge, removed with
    // their faces to kept neighbours re-emitted as wall faces.
    int invalidSliverRemovals = 0;
    // Manifold-or-refuse merges. A merge-component whose
    // WOULD-BE merged cell's full face set (wall + internal) fails the
    // integer-topology manifoldness test (every edge incidence exactly 2,
    // no pinch vertices) is REFUSED -- its sliver members fall back to the
    // same removal path as `mergeFallbackRemovals`, but counted separately
    // (this is a distinct mechanism: a merge that COULD have found an
    // anchor but produced an invalid polyhedron, vs. one that never found
    // any kept neighbour at all).
    int mergeRefusedTopology = 0;
    // Convex-or-refuse merges (the concave-crease family): a merge-component whose WOULD-BE
    // merged cell is manifold but places its centroid OUTSIDE one of its
    // own facets' outward half-space (checkMesh's "incorrectly oriented
    // face pyramid" criterion) is likewise REFUSED -- the members stay
    // their own valid, convex cut cells. Counted separately from
    // mergeRefusedTopology because the mechanism is different (geometry,
    // not topology: the concave-crease folded-facet family). A pure sign
    // test, no epsilon; see Cutter.cpp::isConvexMergedCell.
    int mergeRefusedConvexity = 0;
    // Merges refused because the merged polyhedron failed the
    // checkMesh-mimicking geometry tests (see cellGeometryIsValid).
    int mergeRefusedGeometry = 0;
    // Summed area of the OWN wall-patch faces of every refused-topology
    // sliver (see mergeRefusedTopology): makes visible the trade a refused
    // merge makes -- it KEEPS the sliver as its own small cell rather than
    // growing the solid into it (growing the solid regresses
    // watertightness -- see Cutter.cpp), so this area is wall kept at
    // small-cell resolution, not solid growth.
    double mergeRefusedWallArea = 0.0;
    // Merged-sliver skew gate. A manifold-AND-convex merge candidate can
    // still be geometrically THIN enough that its shared internal faces to
    // outside neighbours score high checkMesh face skewness (the
    // "sliver-column" family). Bound = 3.0
    // (Cutter.cpp mergeSlivers, `kMergeSkewBound`, full provenance and the
    // ctest offset_twospheres_layers collateral-regression story there):
    // a static histogram (on
    // win_fish_coarse_skew2) suggested 3.5 from the naive gap, but that
    // measurably left 13/16 residual skewFaces against real checkMesh (a
    // group's true final skew depends on refusal decisions made for its
    // neighbours too -- an interaction a static histogram cannot see).
    // 3.0 clears the isolated window (0 skewFaces, max skewness 3.93333
    // OK) and is the LOWEST bound that does not regress any existing
    // test: it shrinks the full bm_layers_fishpassage_coarse parent's
    // family from 45 faces/4.72889 to a disclosed partial residual of 7
    // faces/4.04765 -- going lower (2.5) fully clears the parent too but
    // reselects merges on ctest's offset_twospheres_layers that shift the
    // layer-quality gate's stack removal enough to break that
    // test's wall-area sanity band. Honest partial fix, not forced.
    // `mergeReselectedSkew`: a 2-member (single sliver + single anchor)
    // group whose ORIGINAL pairing failed the skew bound found a DIFFERENT
    // admissible partner (another kept, non-sliver, not-already-merged
    // neighbour of the sliver) that passes both convexity and the skew
    // bound -- the merge still happens, just absorbed in a different
    // direction. `mergeRefusedSkew`: a group that failed the skew bound and
    // had no working reselection (either a size>=3 chain -- reselection is
    // not enumerated/attempted there, see Cutter.cpp -- or a
    // 2-member group with no alternative admissible partner at all, or none
    // of its alternatives passed) falls back to the SAME refuse treatment
    // as mergeRefusedConvexity (members stay their own untouched cells).
    int mergeReselectedSkew = 0;
    int mergeRefusedSkew = 0;
    // `mergeKeptOverSkew`: a 2-member group between kMergeSkewBound and
    // kMergeSkewFallback, kept because its sliver would be deleted unmerged.
    int mergeKeptOverSkew = 0;
    // Layer-feasibility detection: populated by the
    // caller from Layers.cpp's LayerStats after applyLayers runs, same
    // "populated by the caller" pattern as residualMean/residualMax.
    double infeasibleWallArea = 0.0;
    // Connectivity drop. Cells removed because they
    // are not reachable, across shared internal faces, from the cell
    // containing `locationInMesh`; discardedComponents counts the
    // distinct unreachable connected components dropped (0 when the
    // kept set was already one component).
    int disconnectedCellsDropped = 0;
    int discardedComponents = 0;
    // Fluid volume in those discarded components. The flow domain the
    // written mesh describes is the STL's minus this -- disclosed with a
    // per-component bbox on stdout and a loud stderr warning, at parity
    // with the layer stage's offsetCutDisconnected* report.
    double discardedVolume = 0.0;
    // Smoothed field: populated by the caller from
    // Layers.cpp's LayerStats, same "populated by the caller" pattern as
    // residualMean/residualMax above. All zero unless smoothRadius
    // is set.
    long smoothGradFallback = 0;
    long smoothMovedFaces = 0;
    double smoothMovedArea = 0.0;
    double smoothResidWeightedSum = 0.0;
    double smoothResidAreaSum = 0.0;
    double smoothResidMaxOverH = 0.0;
    int smoothSealedRegionCount = 0;
    double smoothSealedRegionArea = 0.0;
    // The per-region breakdown behind the two numbers above, largest
    // area first. Empty unless smoothing is on. Count and summed area alone
    // cannot say WHERE the smoothing changed inside/outside, nor whether the
    // same passages seal as the radius grows -- which is the whole question the
    // sealed-region instrument was added to answer.
    std::vector<SealedRegionRecord> smoothSealedRegionDetail;
    // Land-on-phi_r=0: bisection-bracket fallback counts
    // and the achieved stack-thickness distribution's summary,
    // populated by the caller from Layers.cpp's LayerStats, same
    // pattern as the smoothing fields above. All zero/default unless
    // smoothing is active.
    long smoothLandingBisectFallback = 0;
    long smoothLandingDegenerateGradFallback = 0;
    long smoothLandingNoBracketFallback = 0;
    long smoothLandingEvalSum = 0;
    long smoothLandingEvalMax = 0;
    long smoothLandingCandidateCount = 0;
    long stackThicknessCount = 0;
    double stackThicknessMin = 0.0;
    double stackThicknessMean = 0.0;
    double stackThicknessP1 = 0.0;
    long stackThicknessBelowHalfT = 0;
    long stackThicknessBelowQuarterT = 0;
};

// Cuts `base` (a general convex-polyhedron mesh)
// against the shared vertex/edge data in `cd`, producing a compacted
// mesh ready for FoamWriter. Cut and STL-solid-side faces are assigned
// to one of `wallPatchNames` (ordinal == `Triangle::solidId`, dict
// declaration order, appended after `base`'s existing patches): a face
// carrying at least one cut-edge intercept goes to the STL that produced
// the majority of its intercepts (ties -> lowest ordinal); a boundary
// face with no intercepts of its own (its cell became a boundary only
// because a neighbour was fully removed) goes to the STL nearest its
// centroid, via `tris` (the combined triangle soup, `solidId`-tagged).
// Every entry of `wallPatchNames` is emitted as its own patch even if it
// ends up with zero faces (deterministic patch ordering).
// `originCell` (optional): if non-null, filled with one entry
// per output cell giving its originating `base` cell index. Threading
// this cell id through the cutter is what lets
// propagateLevelsThroughCut (Refine.cpp) do an O(1) lookup per output
// cell instead of an O(n_out * n_base) AABB scan -- measured to be
// prohibitively slow at 1M-cell scale.
// `volFracThreshold`: the sliver-drop
// threshold below which a cut cell is snapped fully solid (removed).
// Defaults to `kSmallVolFrac` (1e-3) --
// only the offset-cut path in main.cpp passes a different value
// (`absorbVolFrac`, default 0.3, the design's sliver-sink decision for
// the layer interface). The cutting/classification algorithm itself is
// unaffected; only which volume fraction counts as "negligible" changes.
// `keepSliversForMerge`: when true, a cell that would
// otherwise be dropped by the `volFracThreshold` sliver rule is instead
// KEPT (as a normal cut cell, full geometry) and flagged in
// `sliverFlagsOut` (resized to the output cell count, indexed the same
// as every other per-output-cell array). Default false;
// only `runPipeline`'s offset-cut branch passes true. The actual MERGE
// (deleting shared faces, unioning sliver cells into a kept neighbour)
// is a separate post-pass, `mergeSlivers` below -- cutMesh itself only
// answers "is this cell a sliver", it does not touch topology.
GeneratedMesh cutMesh(const GeneratedMesh& base, const CutData& cd,
                      const std::vector<std::string>& wallPatchNames, const std::vector<Triangle>& tris,
                      CutStats& stats, std::vector<int>* originCell = nullptr, double volFracThreshold = -1.0,
                      bool keepSliversForMerge = false, std::vector<bool>* sliverFlagsOut = nullptr);

// Sliver merge (offset-cut path only): merges every cell
// flagged in `isSliver` (same indexing as `mesh`'s cells, e.g. from
// cutMesh's `sliverFlagsOut`) into a kept (non-sliver) neighbour --
// deleting the shared face(s) between them and unioning the two cells'
// face sets into one polyhedron. Rules:
//  - merge target = the kept neighbour sharing the LARGEST-area face
//    (tie -> lowest cell id);
//  - a connected component of slivers (adjacent via an internal face)
//    agglomerates via union-find into ONE anchor: the kept neighbour of
//    ANY component member sharing the largest-area face (tie -> lowest
//    cell id);
//  - a component with NO kept neighbour at all falls back to REMOVAL
//    (the pre-merge treatment) -- counted in `stats.mergeFallbackRemovals`.
// Mutates `mesh`/`cellLevel`/`pointLevel` in place (cells/points may be
// removed and renumbered, exactly like cutMesh's own compaction).
void mergeSlivers(GeneratedMesh& mesh, std::vector<int>& cellLevel, std::vector<int>& pointLevel,
                   const std::vector<bool>& isSliver, CutStats& stats, const Vec3& locationInMesh);

// Factored out of dropDisconnectedCells so
// debug-only dump sites (NINJA_MERGE_DEBUG_DIR in mergeSlivers,
// NINJA_PRELAYERS_DUMP in main.cpp) can compute, redundantly and at
// whatever mesh state exists at dump time, the SAME "would this cell
// survive dropDisconnectedCells" mask dropDisconnectedCells itself uses --
// without duplicating the seed-location + BFS logic. Locates the cell
// containing `locationInMesh` (same containment test/fallback
// dropDisconnectedCells uses) and flood-fills across shared internal
// faces; returns one entry per CURRENT `mesh.nCells()`, true = reachable
// (would be kept). Throws std::runtime_error under the same "no cell
// contains locationInMesh" condition dropDisconnectedCells does.
std::vector<char> reachableCellsFromLocation(const GeneratedMesh& mesh, const Vec3& locationInMesh);

// Flood-fills the kept cell set across shared
// internal faces starting from the cell that contains `locationInMesh`
// (a half-space containment test against each candidate cell's own
// faces, oriented outward from its owner side; ties/edge cases fall
// back to nearest-centroid) and drops every cell not reached --
// dropped-region faces (internal and boundary) are deleted, kept faces
// are remapped/compacted exactly like `mergeSlivers`' own compaction.
// A no-op (mesh untouched, stats zeroed) when every cell is already
// reachable -- the byte-identical guarantee for the ordinary case.
// Throws std::runtime_error if no kept cell contains `locationInMesh`
// (hard dict/geometry error: fail, do not guess).
// `extraCellFields`: further per-cell int fields (e.g. stackId /
// layerIndex) that must be renumbered by the SAME compaction. A field whose
// length does not match the pre-drop cell count is left untouched.
void dropDisconnectedCells(GeneratedMesh& mesh, std::vector<int>& cellLevel, std::vector<int>& pointLevel,
                            const Vec3& locationInMesh, CutStats& stats,
                            const std::vector<std::vector<int>*>& extraCellFields = {});


// PLANARIZE THE LAYER TOP.
//
// A slightly non-planar cut face is acceptable in general, since OpenFOAM
// tolerates mildly non-planar faces. That is true for a WALL face and
// FALSE for the layer top, which is the seed every prism stack is built on:
// a prism whose base is bent by more than its own first-layer thickness
// cannot be valid, and OpenFOAM's corner-anchored tet decomposition -- the
// exact predicate the layer-quality gate uses -- is what rejects it.
// MEASURED on bm_layers_fishpassage_coarse's pre-layer mesh: 83.67% of the
// 206,269 layer-top faces are already planar to machine precision, 16.33%
// are not, and 0.11% (~227 faces) are bent by more than the first layer
// thickness d1 = 0.00711 -- against 337 condemned stacks.
//
// Fix: split each non-planar boundary face on the named patches into planar
// triangles by fanning from an existing vertex. Properties that make this
// cheap and safe, all deliberate:
//   * NO NEW POINTS. A fan uses existing vertices, so the march's front
//     point set, positions and motion are completely unchanged; only the
//     face list (and hence the front adjacency graph, which gains the fan
//     diagonals) differs.
//   * NO WATERTIGHTNESS RISK. These are BOUNDARY faces with a single owner
//     cell, and a fan adds only interior diagonals -- the face's perimeter
//     loop is untouched, so every neighbouring cell still sees exactly the
//     same shared edges. This is unlike every other change to the cut path,
//     which must satisfy the A+parity shared-data contract.
//   * The owner cell's volume shifts by the warp it was previously
//     misrepresenting, which is a correction: OpenFOAM already evaluates a
//     non-planar face by triangulating it.
// The apex is chosen to maximize the minimum triangle quality among the
// apexes that leave the owner cell valid (see below), and is a pure function
// of the owner cell's own geometry. Non-convex faces are NOT split (a fan
// could escape the polygon) and are counted separately.
//
// `warpTol` is absolute; faces with max deviation from their best-fit plane
// at or below it are left alone. Returns the number of faces split; the
// out-params report triangles added and non-convex faces skipped.
// The apex is admissible only if the OWNER CELL, rebuilt with this fan in
// place of the face, still passes checkMesh's own per-face tests
// (`cellGeometryIsValid`): pyramid orientation about the volume-weighted cell
// centre, a usable vertex-anchored tet decomposition, and relative-area
// non-degeneracy. Quality decides among the admissible apexes.
// `apexPyramidFallback` (optional): faces where NO apex was admissible, so the
// face was LEFT UNSPLIT (warped but valid, which this cutter's founding design
// already permits) rather than split into a fan checkMesh would reject. See
// Cutter.cpp.
// THE CHECKMESH TRANSCRIPTION, and the mesher's core per-cell post-condition.
// `loops` is a cell's COMPLETE emitted face set, every loop wound OUTWARD from
// the cell. True iff every face is non-degenerate (relative area), has a
// positively-oriented pyramid about the volume-weighted cell centre, has an
// all-negative FACE-CENTRE fan (polyMeshTetDecomposition::checkFaceTets' first
// loop -- the only tet test checkMesh applies to a BOUNDARY face, and the one a
// non-convex polygon fails) and admits a usable VERTEX-anchored tet
// decomposition (findBasePoint, checkFaceTets' internal-face half). Transcribed
// from OpenFOAM's sources, not paraphrased; see Cutter.cpp. Exposed for the
// `cell_geometry` unit gate.
bool cellGeometryIsValid(const std::vector<std::vector<int>>& loops, const std::vector<Vec3>& pts);

// Opt-in (`NINJA_GEOM_AUDIT=1`) whole-mesh run of the checkMesh transcription
// `cellGeometryIsValid` at a named pipeline point, reporting the offending
// cells and their failing faces on stderr. No-op and byte-identical when the
// variable is unset; pure diagnostics, used for stage attribution.
void auditCellGeometry(const GeneratedMesh& mesh, const char* label);

int planarizeBoundaryFaces(GeneratedMesh& mesh, const std::vector<std::string>& patchNames,
                           double warpTol, int* trianglesAdded, int* nonConvexSkipped,
                           int* apexPyramidFallback = nullptr);


// COPLANAR-FLAP REPAIR.
//
// Fixes a half-grid defect family measured on a rotated-cube case
// (offset_rotcube geometry):
// 24 wall faces with a NEGATIVE face pyramid (12 at z=0.25, 12 at z=0.75 --
// exactly the two ON planes), which the layer march then seeded stacks on and
// the layer-quality gate condemned, cascading into 420 dropped faces.
//
// What they are. `vertexSolid` is a GLOBAL IN/not-IN flag, so a grid vertex
// lying exactly on some OTHER sheet of the surface reads as fluid (ON) even
// where the sheet cutting THIS cell puts it strictly on the solid side. Every
// edge from such a vertex to an IN vertex then reads as separating, the
// half-grid rule snaps its intercept back onto the vertex itself, and the
// vertex is injected into the cut loop. The loop closes but FOLDS FLAT: the
// sub-arc around the injected vertex is coplanar with, and inside, a grid face
// the cell keeps whole. MEASURED: 48 such faces on rotcube, each exactly three
// of the four corners of an internal face of its own cell.
//
// Why the obvious repairs are wrong, both MEASURED:
//   * Flipping the offending faces opens
//     the same cells -- checkMesh "max cell openness 0.5", precisely the flap's
//     area ratio. The cells are NOT inside-out; all 24 have exactly the right
//     positive volume, and closure held via the flap all along.
//   * Splicing the injected vertex out of the cut loop opens them too (measured
//     openness 0.174 against a 4.5e-17 baseline): the flap is load-bearing,
//     because the host grid face is still present at FULL size on the far side.
//
// What this does instead. The flap and the host face are the same surface
// counted twice, on the wrong sides. Lift the flap arc off the boundary face
// and RE-PARENT it to the cell on the far side of the host, deleting the arc's
// interior vertices from both loops. Both partitions are then exact --
// boundary = remainder + flap, host = hostLoop + flap -- so every affected
// cell's face-area sum is algebraically unchanged and BOTH stay watertight,
// while the flap now bounds the cell it actually bounds and its pyramid is
// positive. Volumes are unchanged. Only faces already violating the pyramid
// invariant are touched, so a healthy mesh is left byte-identical.
//
// `patchNames` limits the search to the layered wall patches. Returns the
// number of flaps re-parented.
// Counts cells whose emitted face set is not a closed polyhedron: fewer
// than four faces, or an undirected edge of the cell's boundary shared by
// other than exactly two of its faces. This is the in-mesher form of the
// invariant `Benchmarks/check_integrity.py` measures on the written mesh.
//
// Zero is the contract for the core mesh. cutMesh enforces it on the cells
// it emits (see its closed-polyhedron post-condition); this function is
// what lets the stages that RE-WRITE that mesh afterwards -- mergeSlivers,
// repairCoplanarFlaps, planarizeBoundaryFaces -- be held to the same
// contract instead of being trusted. Reported per stage under
// NINJA_CORE_INTEGRITY.
int countNonClosedCells(const GeneratedMesh& mesh);

int repairCoplanarFlaps(GeneratedMesh& mesh, const std::vector<std::string>& patchNames);

// FACE-PYRAMID TRIPWIRE (diagnostic; never mutates the mesh).
//
// Counts face/cell incidences violating the invariant OpenFOAM's face-pyramid
// check tests -- dot(Sf, faceCentre - cellCentre) must be strictly positive on
// the owner side and strictly negative on the neighbour side. The return value
// is directly comparable to checkMesh's "Error in face pyramids: N faces are
// incorrectly oriented".
//
// Deliberately NOT a repair. Flipping the offending faces here instead is
// MEASURED to be wrong: on
// the rotated-cube half-grid case above the 24 offenders are
// coplanar FLAPS -- each is exactly three of the four corners of a grid face
// its own cell keeps whole, so it duplicates half of an existing face at zero
// volume. The cells are NOT inside-out (all 24 have exactly the right positive
// volume, 2.0833e-05); closure holds only by cancellation between the flap and
// its host face, which is why flipping the flap opened the cell (openness 0.5,
// precisely the area ratio). The defect is in the cut LOOP, and cutMesh's
// coplanar-flap pass is where it is now fixed. This function exists so that the
// next member of the family is reported rather than silently patched over.
int countBadFacePyramids(const GeneratedMesh& mesh);

} // namespace ninja
