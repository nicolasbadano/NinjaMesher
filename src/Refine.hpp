// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "BaseMesh.hpp"
#include "CutData.hpp"
#include "Stl.hpp"

namespace ninja {

// One `refinement { zoneN { ... } }` region, including
// "surface" distance-mode regions.
// A leaf is refined toward `level` iff its centroid lies inside the
// region (box: componentwise min<=c<=max; sphere: |c-centre|<=radius;
// surface: distance to the named STL's triangles <= `distance`).
struct RefineRegion {
    std::string type; // "box", "sphere", or "surface"
    Vec3 min;
    Vec3 max;
    Vec3 centre;
    double radius = 0.0;
    int level = 0;
    // surface-mode only: `stlKey` is the geometry dict's raw filename key
    // (e.g. "box.stl") -- resolved by the caller into a triangle list via
    // `refineLocal`'s `stlTriangles` map. `distance` is the band radius.
    std::string stlKey;
    double distance = 0.0;
};

// --- MPI: rank-local refinement + rank-0 assembly ------------------------
//
// Proper domain decomposition: each rank
// owns ONLY its static 3D base-cell block (parDims3) plus a
// one-base-cell halo. `refineLocal` runs the mark -> split -> grade
// loop on the rank-local octree, exchanging halo LEAF LEVELS with
// face-neighbour ranks each grading iteration (three sequential
// axis phases whose slabs include previously-received halo extents, so
// edge/corner-diagonal halo cells arrive through face-neighbour
// messages only) and terminating on an Allreduce fixpoint. The
// conformance-splice corner keys are DERIVED locally from the exchanged
// halo levels (a quad's corners are a pure function of leaf geometry),
// which realizes the corner-key exchange with the same
// one-round information flow and no separate message type.
//
// `refineLocal` also enumerates, per rank, the block-local cut
// entities: the points and (uniquely-owned) edges of this rank's owned
// faces' spliced loops, exactly the sub-segments the serial mesh's
// `collectEdges` would produce there -- so classification/intercepts
// run rank-local and the gathered CutData is bit-identical to the
// serial one (every value is a pure function of global position).
//
// `assembleRefined` is rank-0-only: it rebuilds the global final tree
// from the gathered leaves and runs the UNCHANGED serial assembly
// (leaf numbering, quad emission, conformance splice, face adjacency,
// point ids) -- a pure function of the final leaf set, so its output is
// byte-identical for every rank count.

// One final octree leaf in GLOBAL fine-grid units (POD wire record for
// the leaf gather; `level` is the leaf's depth, extent R>>level).
struct LeafRec {
    int i0 = 0, j0 = 0, k0 = 0;
    int level = 0;
};

// Global fine-grid LATTICE flat index (point key), shared by
// refineLocal / assembleRefined / main's gather reconstruction:
// i + j*(nxF+1) + k*(nxF+1)*(nyF+1), with nxF = nx*R etc.
inline long long fineLatticeFlat(int i, int j, int k, int nxF, int nyF) {
    return static_cast<long long>(i) +
           static_cast<long long>(j) * (nxF + 1) +
           static_cast<long long>(k) * (static_cast<long long>(nxF + 1) * (nyF + 1));
}

// Rank-local refinement + cut-entity enumeration result.
struct LocalRefine {
    int R = 1; // 2^maxLevel
    // This rank's OWNED leaves of the final (globally graded) tree.
    std::vector<LeafRec> ownedLeaves;
    // Block-local cut entities (empty when needCutEntities=false):
    std::vector<Vec3> points;           // classification points (coords)
    std::vector<long long> pointFlat;   // global lattice flat per point
    // Level of the FINEST leaf touching each point (the 8 fine cells
    // around its lattice position). Read off the final graded tree
    // including the halo, so it is a pure function of global position --
    // identical on every rank that enumerates the point.
    std::vector<int> pointFinestLevel;
    std::vector<EdgeKey> edges;         // OWNED edges (local point-index
                                        // pairs; each global edge is
                                        // owned by exactly one rank)
};

// Collective (all ranks). `regions` empty runs the same code with R=1
// and no marking work; point coordinates then use generateBaseMesh's
// exact expression (min + (max-min)*(i/n)) instead of the refined
// path's (min + (max-min)*i/(n*R)) so gathered values stay
// bit-identical to whichever serial mesh rank 0 assembles.
LocalRefine refineLocal(const MeshConfig& cfg, const std::vector<RefineRegion>& regions,
                        const std::unordered_map<std::string, std::vector<Triangle>>& stlTriangles,
                        bool needCutEntities);

// Rank-0-only serial assembly from the gathered global leaf set.
// `pointFlat[p]` gives mesh point p's global lattice flat (for mapping
// the gathered classification/intercept records onto mesh point ids).
struct AssembledRefine {
    GeneratedMesh mesh;
    std::vector<int> cellLevel;
    std::vector<int> pointLevel;
    std::vector<long long> pointFlat;
};
AssembledRefine assembleRefined(const MeshConfig& cfg, const std::vector<LeafRec>& leaves, int R);

// Per-level cell counts and the max face-adjacent cellLevel difference
// (2:1 grading verification), computed from a mesh's internal-face
// adjacency + its cellLevel array. Used by --refine-stats.
struct RefineStats {
    int maxAdjacentLevelDiff = 0;
    std::vector<int> countByLevel; // index = level
};
RefineStats computeRefineStats(const GeneratedMesh& mesh, const std::vector<int>& cellLevel);

// Propagates cellLevel/pointLevel from the pre-cut refined mesh through
// the (unmodified) cutter's output. Cell identity survives the cut
// (cutMesh only drops/keeps cells, never merges or splits one).
// `originCell` (optional -- see Cutter.hpp's `cutMesh`): when
// provided (one entry per output cell, its pre-cut cell index), each
// kept output cell's level is an O(1) lookup via that id. When null
// (back-compat), falls back to the original O(n_out * n_precut) AABB
// containment scan -- MEASURED to blow past the 10-minute gate at
// bm_scale_refined's 1M-cell scale, which is why every call site should
// pass `originCell` now. Point identity is exact-match: pre-cut points
// keep their level (looked up by exact coordinate); a cut-intercept
// point (not present pre-cut) takes max(level of the edge's two
// endpoints), per `cd`.
struct PostCutLevels {
    std::vector<int> cellLevel;
    std::vector<int> pointLevel;
};
PostCutLevels propagateLevelsThroughCut(const GeneratedMesh& preCut,
                                         const std::vector<int>& preCutCellLevel,
                                         const std::vector<int>& preCutPointLevel, const CutData& cd,
                                         const GeneratedMesh& cut, const std::vector<int>* originCell = nullptr);

} // namespace ninja
