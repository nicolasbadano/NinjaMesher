// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#include "BaseMesh.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <unordered_map>

namespace ninja {

MergedPatches mergePatchSpecs(const std::vector<PatchSpec>& specs) {
    MergedPatches out;
    std::unordered_map<std::string, std::size_t> nameToOrdinal;
    for (const PatchSpec& spec : specs) {
        auto it = nameToOrdinal.find(spec.name);
        if (it == nameToOrdinal.end()) {
            const std::size_t ord = out.patches.size();
            nameToOrdinal.emplace(spec.name, ord);
            PatchInfo info;
            info.name = spec.name;
            info.type = spec.type;
            out.patches.push_back(info);
            out.sideKeyToOrdinal[spec.sideKey] = ord;
        } else {
            out.sideKeyToOrdinal[spec.sideKey] = it->second;
            if (out.patches[it->second].type != spec.type) {
                std::cerr << "warning: domain side '" << spec.sideKey << "' declares type '" << spec.type
                          << "' but patch '" << spec.name << "' was already declared as '"
                          << out.patches[it->second].type << "'; keeping the first\n";
            }
        }
    }
    return out;
}


namespace {

// Hash for a sorted small vector of point indices ("sorted-point-key"),
// used to match faces shared between adjacent cells. This is the same
// mechanism the future STL cutter will use for cross-cell face
// matching, exercised here on the trivial structured-hex case.
struct KeyHash {
    std::size_t operator()(const std::vector<int>& key) const noexcept {
        std::size_t h = key.size();
        for (int v : key) {
            h ^= static_cast<std::size_t>(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return h;
    }
};

std::vector<int> sortedKey(const std::array<int, 4>& pts) {
    std::vector<int> k(pts.begin(), pts.end());
    std::sort(k.begin(), k.end());
    return k;
}

int pointIndex(int i, int j, int k, int nx, int ny) {
    return i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
}

int cellIndex(int i, int j, int k, int nx, int ny) {
    return i + j * nx + k * nx * ny;
}

// One of the 6 local faces of a hex cell, described generically so the
// same loop drives both adjacency-building passes.
struct LocalFace {
    std::array<int, 4> points; // outward-from-this-cell order
    const char* sideKey;       // e.g. "xmin"; nullptr if this local face
                                // cannot be a domain boundary for this cell
};

std::array<LocalFace, 6> localFaces(int i, int j, int k, int nx, int ny, int nz) {
    const int p000 = pointIndex(i, j, k, nx, ny);
    const int p100 = pointIndex(i + 1, j, k, nx, ny);
    const int p010 = pointIndex(i, j + 1, k, nx, ny);
    const int p110 = pointIndex(i + 1, j + 1, k, nx, ny);
    const int p001 = pointIndex(i, j, k + 1, nx, ny);
    const int p101 = pointIndex(i + 1, j, k + 1, nx, ny);
    const int p011 = pointIndex(i, j + 1, k + 1, nx, ny);
    const int p111 = pointIndex(i + 1, j + 1, k + 1, nx, ny);

    std::array<LocalFace, 6> faces;
    faces[0] = LocalFace{{p000, p001, p011, p010}, i == 0 ? "xmin" : nullptr};
    faces[1] = LocalFace{{p100, p110, p111, p101}, i == nx - 1 ? "xmax" : nullptr};
    faces[2] = LocalFace{{p000, p100, p101, p001}, j == 0 ? "ymin" : nullptr};
    faces[3] = LocalFace{{p010, p011, p111, p110}, j == ny - 1 ? "ymax" : nullptr};
    faces[4] = LocalFace{{p000, p010, p110, p100}, k == 0 ? "zmin" : nullptr};
    faces[5] = LocalFace{{p001, p101, p111, p011}, k == nz - 1 ? "zmax" : nullptr};
    return faces;
}

// Adjacency record built in pass 1, keyed by sorted-point-key.
struct AdjRecord {
    int firstCell = -1;
    int secondCell = -1; // -1 until a second cell registers this face
};

} // namespace

void buildCellFaces(GeneratedMesh& mesh, int nCells) {
    std::vector<int> count(static_cast<std::size_t>(nCells), 0);
    const int nF = mesh.faces.size();
    for (int f = 0; f < nF; ++f) {
        ++count[static_cast<std::size_t>(mesh.faces.owner[static_cast<std::size_t>(f)])];
        const int nb = mesh.faces.neighbour[static_cast<std::size_t>(f)];
        if (nb != -1) {
            ++count[static_cast<std::size_t>(nb)];
        }
    }
    mesh.cellFaceOffsets.assign(static_cast<std::size_t>(nCells) + 1, 0);
    for (int c = 0; c < nCells; ++c) {
        mesh.cellFaceOffsets[static_cast<std::size_t>(c) + 1] =
            mesh.cellFaceOffsets[static_cast<std::size_t>(c)] + count[static_cast<std::size_t>(c)];
    }
    mesh.cellFaces.assign(static_cast<std::size_t>(mesh.cellFaceOffsets[static_cast<std::size_t>(nCells)]), -1);
    std::vector<int> cursor(mesh.cellFaceOffsets.begin(), mesh.cellFaceOffsets.end() - 1);
    for (int f = 0; f < nF; ++f) {
        const int own = mesh.faces.owner[static_cast<std::size_t>(f)];
        mesh.cellFaces[static_cast<std::size_t>(cursor[static_cast<std::size_t>(own)]++)] = f;
        const int nb = mesh.faces.neighbour[static_cast<std::size_t>(f)];
        if (nb != -1) {
            mesh.cellFaces[static_cast<std::size_t>(cursor[static_cast<std::size_t>(nb)]++)] = f;
        }
    }
}

GeneratedMesh generateBaseMesh(const MeshConfig& cfg) {
    const int nx = cfg.nx;
    const int ny = cfg.ny;
    const int nz = cfg.nz;
    if (nx <= 0 || ny <= 0 || nz <= 0) {
        throw std::runtime_error("BaseMesh error: n must be positive in all directions");
    }

    GeneratedMesh mesh;

    // Points.
    mesh.points.resize(static_cast<std::size_t>(nx + 1) * (ny + 1) * (nz + 1));
    for (int k = 0; k <= nz; ++k) {
        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
                const double x = cfg.min.x + (cfg.max.x - cfg.min.x) * (static_cast<double>(i) / nx);
                const double y = cfg.min.y + (cfg.max.y - cfg.min.y) * (static_cast<double>(j) / ny);
                const double z = cfg.min.z + (cfg.max.z - cfg.min.z) * (static_cast<double>(k) / nz);
                mesh.points[pointIndex(i, j, k, nx, ny)] = Vec3{x, y, z};
            }
        }
    }

    const int nCells = nx * ny * nz;

    // Pass 1: build the face->cells adjacency (FaceCellAdjacency) via
    // sorted-point-key matching, during mesh generation (not a writer
    // post-pass).
    std::unordered_map<std::vector<int>, AdjRecord, KeyHash> adjacency;
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const int c = cellIndex(i, j, k, nx, ny);
                for (const LocalFace& lf : localFaces(i, j, k, nx, ny, nz)) {
                    const std::vector<int> key = sortedKey(lf.points);
                    AdjRecord& rec = adjacency[key];
                    if (rec.firstCell == -1) {
                        rec.firstCell = c;
                    } else {
                        rec.secondCell = c;
                    }
                }
            }
        }
    }

    // Group boundary patches in dict declaration order, by patch NAME
    // (see mergePatchSpecs).
    const MergedPatches merged = mergePatchSpecs(cfg.patches);
    const std::vector<PatchInfo>& mergedPatches = merged.patches;
    const std::unordered_map<std::string, std::size_t>& sideKeyToPatchOrdinal = merged.sideKeyToOrdinal;
    std::vector<FaceStore> boundaryByPatch(mergedPatches.size());

    // Pass 2: emit faces exactly once, from the cell that owns them
    // (lower cell index — or the sole cell, for boundary faces), in
    // ascending cell order. Internal faces collected separately so
    // they can be sorted into upper-triangular order.
    FaceStore internalFaces;
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const int c = cellIndex(i, j, k, nx, ny);
                for (const LocalFace& lf : localFaces(i, j, k, nx, ny, nz)) {
                    const std::vector<int> key = sortedKey(lf.points);
                    const AdjRecord& rec = adjacency.at(key);
                    const bool isBoundary = (rec.secondCell == -1);
                    const int owner = isBoundary ? rec.firstCell
                                                  : std::min(rec.firstCell, rec.secondCell);
                    const int neighbour = isBoundary ? -1
                                                      : std::max(rec.firstCell, rec.secondCell);
                    if (owner != c) {
                        continue; // emitted when the owner cell is visited
                    }

                    if (isBoundary) {
                        if (lf.sideKey == nullptr) {
                            throw std::runtime_error(
                                "BaseMesh internal error: unmatched face is not on a domain boundary");
                        }
                        auto it = sideKeyToPatchOrdinal.find(lf.sideKey);
                        if (it == sideKeyToPatchOrdinal.end()) {
                            throw std::runtime_error(
                                std::string("BaseMesh error: no patch declared for domain side '") +
                                lf.sideKey + "'");
                        }
                        boundaryByPatch[it->second].append(lf.points.data(), 4, owner, neighbour,
                                                            static_cast<int>(it->second));
                    } else {
                        internalFaces.append(lf.points.data(), 4, owner, neighbour, -1);
                    }
                }
            }
        }
    }

    // Internal faces are already emitted in owner-ascending order, and
    // for a fixed owner the local-face traversal order (-x,+x,-y,+y,-z,+z)
    // yields ascending neighbour indices for the structured grid.
    // Sort explicitly anyway to make the upper-triangular guarantee
    // independent of traversal order.
    internalFaces.stableSortByOwnerNeighbour();

    mesh.nInternalFaces = internalFaces.size();
    mesh.faces = std::move(internalFaces);

    for (std::size_t p = 0; p < mergedPatches.size(); ++p) {
        PatchInfo info = mergedPatches[p];
        info.startFace = mesh.faces.size();
        info.nFaces = boundaryByPatch[p].size();
        mesh.faces.appendAll(boundaryByPatch[p]);
        mesh.patches.push_back(info);
    }

    // Wire cells' face adjacency (polyhedral-ready: a cell is a list of
    // face indices, not a fixed 8-point hex record).
    buildCellFaces(mesh, nCells);

    return mesh;
}

} // namespace ninja
