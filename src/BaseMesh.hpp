// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace ninja {

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

// One entry of the `patches` dict block, in declaration order.
// sideKey is one of "xmin", "xmax", "ymin", "ymax", "zmin", "zmax".
struct PatchSpec {
    std::string sideKey;
    std::string type; // emitted verbatim into boundary file
    std::string name;
};

struct MeshConfig {
    Vec3 min;
    Vec3 max;
    int nx = 0;
    int ny = 0;
    int nz = 0;
    std::vector<PatchSpec> patches;
};

// Checked narrowing cast for CSR offsets/counts.
// int32 indices are the project convention; totals beyond 2^31-1 are a
// pathological future case and must THROW, never wrap.
inline int checkedInt(std::size_t v) {
    if (v > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("Mesh size error: index total exceeds int32 range");
    }
    return static_cast<int>(v);
}

// Read-only view of one face's point loop (or one cell's face list)
// inside a flat CSR array -- supports range-for and indexing.
struct IntSpan {
    const int* b = nullptr;
    const int* e = nullptr;
    const int* begin() const { return b; }
    const int* end() const { return e; }
    int size() const { return static_cast<int>(e - b); }
    bool empty() const { return b == e; }
    int operator[](int i) const { return b[i]; }
};

// Compressed-sparse-row face list -- flat arrays instead of a
// `std::vector<Face>` of `Face{vector<int> points}`, which costs 2 heap
// blocks + ~48 B header per face and millions of small mallocs at
// real-CAD scale. Faces are appended in
// emission order; `offsets` has size()+1 entries, `points` is the flat
// concatenation of every face's loop, and owner/neighbour/patchId are
// parallel per-face arrays (`patchId` -1 for internal faces, else the
// ordinal into the owning mesh's `patches` list). Robustness rules
// (plan): plain int32 everywhere with checkedInt at the append site, no
// packing tricks, points stay double elsewhere.
struct FaceStore {
    std::vector<int> offsets{0};
    std::vector<int> points; // flat loops
    std::vector<int> owner;
    std::vector<int> neighbour;
    std::vector<int> patchId;

    int size() const { return static_cast<int>(owner.size()); }
    int pointCount(int f) const {
        return offsets[static_cast<std::size_t>(f) + 1] - offsets[static_cast<std::size_t>(f)];
    }
    IntSpan pointsOf(int f) const {
        const int* base = points.data();
        return IntSpan{base + offsets[static_cast<std::size_t>(f)],
                       base + offsets[static_cast<std::size_t>(f) + 1]};
    }
    int* mutablePointsOf(int f) { return points.data() + offsets[static_cast<std::size_t>(f)]; }

    void append(const int* pts, int n, int own, int nei, int pid) {
        points.insert(points.end(), pts, pts + n);
        offsets.push_back(checkedInt(points.size()));
        owner.push_back(own);
        neighbour.push_back(nei);
        patchId.push_back(pid);
    }
    void append(const std::vector<int>& pts, int own, int nei, int pid) {
        append(pts.data(), static_cast<int>(pts.size()), own, nei, pid);
    }
    // Copy face `f` of `src` verbatim.
    void appendFrom(const FaceStore& src, int f) {
        append(src.points.data() + src.offsets[static_cast<std::size_t>(f)], src.pointCount(f),
               src.owner[static_cast<std::size_t>(f)], src.neighbour[static_cast<std::size_t>(f)],
               src.patchId[static_cast<std::size_t>(f)]);
    }
    // Copy face `f` of `src` with overridden owner/neighbour/patchId.
    void appendFrom(const FaceStore& src, int f, int own, int nei, int pid) {
        append(src.points.data() + src.offsets[static_cast<std::size_t>(f)], src.pointCount(f), own, nei, pid);
    }
    // Append every face of `src` in order.
    void appendAll(const FaceStore& src) {
        const int base = checkedInt(points.size());
        points.insert(points.end(), src.points.begin(), src.points.end());
        for (std::size_t f = 1; f < src.offsets.size(); ++f) {
            offsets.push_back(base + src.offsets[f]);
        }
        owner.insert(owner.end(), src.owner.begin(), src.owner.end());
        neighbour.insert(neighbour.end(), src.neighbour.begin(), src.neighbour.end());
        patchId.insert(patchId.end(), src.patchId.begin(), src.patchId.end());
    }
    void clear() {
        offsets.assign(1, 0);
        points.clear();
        owner.clear();
        neighbour.clear();
        patchId.clear();
    }

    // Stable sort by (owner, neighbour) ascending -- the project's
    // "upper-triangular" internal-face order. A stable sort of an index
    // permutation with the identical comparator yields the identical
    // final sequence the former stable_sort over `vector<Face>` did
    // (stability resolves ties by original position in both forms).
    void stableSortByOwnerNeighbour() {
        const int n = size();
        std::vector<int> idx(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) idx[static_cast<std::size_t>(i)] = i;
        std::stable_sort(idx.begin(), idx.end(), [this](int a, int b) {
            if (owner[static_cast<std::size_t>(a)] != owner[static_cast<std::size_t>(b)]) {
                return owner[static_cast<std::size_t>(a)] < owner[static_cast<std::size_t>(b)];
            }
            return neighbour[static_cast<std::size_t>(a)] < neighbour[static_cast<std::size_t>(b)];
        });
        FaceStore sorted;
        sorted.points.reserve(points.size());
        sorted.offsets.reserve(offsets.size());
        sorted.owner.reserve(owner.size());
        sorted.neighbour.reserve(neighbour.size());
        sorted.patchId.reserve(patchId.size());
        for (int i : idx) {
            sorted.appendFrom(*this, i);
        }
        *this = std::move(sorted);
    }
};

struct PatchInfo {
    std::string name;
    std::string type;
    int nFaces = 0;
    int startFace = 0;
};

// Domain sides may legitimately share one patch NAME (the usual
// `ymin`/`ymax`/`zmin` all named `sides`, mirroring how the same model
// is written for blockMesh/snappy). Emitting one boundary patch per
// SIDE in that case repeats the name in constant/polyMesh/boundary,
// which OpenFOAM rejects ("Duplicate boundary patch ... Boundary
// definition is in error"). This collapses the six side specs onto one
// patch per distinct name, in first-declaration order.
struct MergedPatches {
    std::vector<PatchInfo> patches;                               // one per distinct name
    std::unordered_map<std::string, std::size_t> sideKeyToOrdinal; // sideKey -> index into `patches`
};
MergedPatches mergePatchSpecs(const std::vector<PatchSpec>& specs);


// The face->cells adjacency mesh in CSR (flat-array) form -- a per-object
// `vector<Face>`/`vector<Cell>` representation measured ~3.5 KB/cell peak
// on real CAD, versus a few hundred bytes of information content.
// Face order:
// all internal faces first (sorted by owner, ties by neighbour, both
// ascending), then boundary faces grouped by patch in dict declaration
// order. Cell->face adjacency is CSR as well, built by
// buildCellFaces() in ascending face order (owner entry before
// neighbour entry per face -- the exact push order of the former
// per-cell vectors).
struct GeneratedMesh {
    std::vector<Vec3> points;
    FaceStore faces;
    std::vector<int> cellFaceOffsets{0};
    std::vector<int> cellFaces;
    std::vector<PatchInfo> patches;
    int nInternalFaces = 0;

    int nFaces() const { return faces.size(); }
    int nCells() const { return static_cast<int>(cellFaceOffsets.size()) - 1; }
    IntSpan cellFacesOf(int c) const {
        const int* base = cellFaces.data();
        return IntSpan{base + cellFaceOffsets[static_cast<std::size_t>(c)],
                       base + cellFaceOffsets[static_cast<std::size_t>(c) + 1]};
    }
};

// Rebuilds `mesh.cellFaceOffsets`/`cellFaces` for `nCells` cells from
// the faces' owner/neighbour arrays: two passes (count, fill), filling
// in ascending face-index order with each face's owner entry before its
// neighbour entry -- byte-for-byte the order the former per-cell
// `faceIndices.push_back` loops produced.
void buildCellFaces(GeneratedMesh& mesh, int nCells);

// Generates a structured nx*ny*nz hexahedral base mesh over the box
// [min, max]. Point index convention: i + j*(nx+1) + k*(nx+1)*(ny+1)
// (x fastest). Cell index convention: i + j*nx + k*nx*ny.
GeneratedMesh generateBaseMesh(const MeshConfig& cfg);

} // namespace ninja
