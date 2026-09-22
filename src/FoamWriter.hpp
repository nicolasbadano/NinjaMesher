// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <vector>

#include "BaseMesh.hpp"

namespace ninja {

// Emits an OpenFOAM v2406 ASCII polyMesh (points, faces, owner,
// neighbour, boundary) into <caseDir>/constant/polyMesh/, consuming
// the already-built face->cells adjacency in a single emit pass.
// Throws std::runtime_error on I/O failure.
void writePolyMesh(const GeneratedMesh& mesh, const std::string& caseDir);

// Emits the AMR-seeding `cellLevel`/`pointLevel` ASCII labelList files
// into <caseDir>/constant/polyMesh/.
// Always called — unrefined meshes get all-zero fields. Must run after
// writePolyMesh (relies on the same directory existing).
void writeLevelFields(const std::vector<int>& cellLevel, const std::vector<int>& pointLevel,
                       const std::string& caseDir);

// Debug tagging: one arbitrary-named ASCII labelList field
// into <caseDir>/constant/polyMesh/<fieldName>, same format as
// cellLevel/pointLevel above (writeLevelFields), so ParaView's
// OpenFOAM reader picks it up the same way. Debug-only use: e.g.
// NINJA_PRELAYERS_DUMP writes `keptFinal` (1/0) here so the
// pre-layers/pre-drop mesh dump can be threshold-filtered to the
// finally-kept component. Must run after writePolyMesh.
void writeCellIntField(const std::vector<int>& values, const std::string& caseDir, const std::string& fieldName);

// Opt-in (NINJA_DIAG_FIELDS) visualization aid: emits four ASCII
// volScalarFields into <caseDir>/0/ so the mesh can be colored in
// ParaView/OpenFOAM without any solver run:
//   refinementLevel -- cellLevel as a scalar
//   isLayer         -- 1 for layer-prism cells, 0 for core cells
//   layerNumber     -- march step that created the cell (-1 for core)
//   cellSize        -- cbrt(cell volume), volumes from the divergence
//                      theorem over the face loops (centroid fan)
// `cellLayerIndex` may be empty (no layers configured): isLayer is
// all-zero and layerNumber all -1 then. Must run after writePolyMesh.
void writeDiagnosticFields(const GeneratedMesh& mesh, const std::vector<int>& cellLevel,
                           const std::vector<int>& cellLayerIndex, const std::string& caseDir);

} // namespace ninja
