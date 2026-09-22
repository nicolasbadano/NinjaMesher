// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "BaseMesh.hpp"

namespace ninja {

// Shared legacy-VTK ASCII polydata writer, used by Layers.cpp and
// NINJA_MERGE_DEBUG_DIR's dumps in Cutter.cpp. Deterministic (no
// pointers, no timestamps); files overwrite on rerun.
//
// A single instance emits AT MOST ONE of POLYGONS / LINES / VERTICES;
// `cellScalars` binds to whichever
// one is populated (checked in that priority order) -- exactly mirrors the
// original Layers.cpp writer when only `polys`/`cellScalars` are used.
struct DebugVtk {
    std::vector<Vec3> pts;

    // Polygon cells (faces): one entry per polygon, values are indices
    // into `pts`.
    std::vector<std::vector<int>> polys;

    // Line cells (edges): one (a,b) pair per line, indices into `pts`.
    std::vector<std::pair<int, int>> lines;

    // Vertex cells (points-as-cells): indices into `pts`.
    std::vector<int> vertexPts;

    // CellData: bound to whichever of polys/lines/vertexPts is non-empty
    // (checked in that order -- a given DebugVtk instance populates
    // exactly one of the three).
    std::vector<std::pair<std::string, std::vector<double>>> cellScalars;

    std::vector<std::pair<std::string, std::vector<double>>> pointScalars;
    std::vector<std::pair<std::string, std::vector<Vec3>>> pointVectors;

    void write(const std::string& path) const;
};

} // namespace ninja
