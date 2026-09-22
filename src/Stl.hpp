// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <vector>

#include "BaseMesh.hpp"

namespace ninja {

struct Triangle {
    Vec3 v0;
    Vec3 v1;
    Vec3 v2;
    // Ordinal of the STL/solid this triangle came from (multi-STL
    // support). readStl() does not set this — callers combining
    // multiple STLs into one triangle soup stamp it per-file after
    // reading (default 0, so single-STL callers are unaffected).
    int solidId = 0;
};

// Reads a binary or ASCII STL file (format detected from the header /
// declared triangle count, no library). No geometric validation is
// performed (input is assumed watertight,
// non-self-intersecting, non-degenerate input).
std::vector<Triangle> readStl(const std::string& path);

} // namespace ninja
