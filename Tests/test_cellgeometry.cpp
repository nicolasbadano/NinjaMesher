// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

// Unit gate on `cellGeometryIsValid`, the mesher's transcription of the
// per-face geometry tests checkMesh applies (Cutter.cpp).
//
// WHY THIS EXISTS, and why it is a UNIT test rather than a window case.
// bm_layers_wfp shipped ONE residual "Error in face tets" face -- a 7-gon
// `fixedWalls` cut face at (3.61 5.19 5.46) -- that the pipeline's own
// post-condition could not see: a whole-mesh audit (NINJA_GEOM_AUDIT)
// reported badCells = 0 after mergeSlivers, after the coplanar-flap repair
// and after planarize, on the very mesh checkMesh then rejected. The
// transcription was incomplete, not the mesh checks mis-wired:
// polyMeshTetDecomposition::checkFaceTets runs a FACE-CENTRE-anchored fan
// over EVERY face and only adds the vertex-anchored findSharedBasePoint test
// for INTERNAL faces, so on a boundary face the old code tested the one thing
// checkMesh does not and skipped the only thing it does.
//
// The defect is a property of ONE cell's face loops, not of the domain, so it
// needs neither the 1.1 M-cell parent nor a grid-aligned sub-box: the two
// polyhedra below are written out by hand and the predicate is called
// directly. That makes the gate exact, sub-millisecond, and immune to the
// parent case being stripped.

#include "Cutter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using ninja::Vec3;

namespace {

int failures = 0;

void expect(bool cond, const char* what) {
    if (!cond) {
        std::printf("FAIL: %s\n", what);
        ++failures;
    } else {
        std::printf("ok:   %s\n", what);
    }
}

// A unit cube, points ordered so that every loop below is wound OUTWARD.
//   0 (0,0,0) 1 (1,0,0) 2 (1,1,0) 3 (0,1,0)   z = 0
//   4 (0,0,1) 5 (1,0,1) 6 (1,1,1) 7 (0,1,1)   z = 1
std::vector<Vec3> cubePoints() {
    return {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
}

std::vector<std::vector<int>> cubeLoops() {
    return {
        {0, 3, 2, 1}, // z = 0, outward = -z
        {4, 5, 6, 7}, // z = 1, outward = +z
        {0, 1, 5, 4}, // y = 0
        {2, 3, 7, 6}, // y = 1
        {0, 4, 7, 3}, // x = 0
        {1, 2, 6, 5}, // x = 1
    };
}

} // namespace

int main() {
    // 1. The plain hex is valid. Guards against a transcription so strict it
    //    rejects the mesh's own bread-and-butter cell.
    expect(ninja::cellGeometryIsValid(cubeLoops(), cubePoints()), "unit cube is valid");

    // 2. THE REGRESSION. Replace the z = 1 quad by a NON-CONVEX hexagon in the
    //    same plane, deep enough that the polygon's area-weighted centroid
    //    falls OUTSIDE the notch edge. Every vertex-anchored fan can still be
    //    fine, so the old (findBasePoint-only) transcription accepted this;
    //    checkMesh's face-centre fan does not, and reports the face under
    //    "Error in face tets". The shape is the same kind of reflex wall facet
    //    the cut and the sliver merge can hand a cell, and the one planarize
    //    declines to split (a vertex fan out of a reflex vertex could leave
    //    the polygon).
    {
        std::vector<Vec3> pts = cubePoints();
        // z = 1 face becomes 4 (0,0,1) 5 (1,0,1) 8 (1,0.9,1) 9 (0.5,0.1,1)
        //                    10 (0,0.9,1) -- a deep notch biting in from +y.
        pts.push_back({1.0, 0.9, 1.0});  // 8
        pts.push_back({0.5, 0.1, 1.0});  // 9  reflex
        pts.push_back({0.0, 0.9, 1.0});  // 10
        std::vector<std::vector<int>> loops = cubeLoops();
        loops[1] = {4, 5, 8, 9, 10};
        // The side walls are left as the cube's: this polyhedron is not
        // closed in the watertight sense, and does not need to be -- the
        // predicate is a per-face test about one cell's loops and its
        // volume-weighted centre, which is exactly what is under test here.
        expect(!ninja::cellGeometryIsValid(loops, pts),
               "cell with a reflex boundary face whose centroid falls outside an edge is REJECTED");
    }

    // 3. The convex counterpart of the same face -- the notch pushed out so the
    //    polygon is convex -- must still be accepted, so the new test is not a
    //    blanket ban on non-quad wall facets.
    {
        std::vector<Vec3> pts = cubePoints();
        pts.push_back({1.0, 0.9, 1.0});  // 8
        pts.push_back({0.5, 1.0, 1.0});  // 9  convex
        pts.push_back({0.0, 0.9, 1.0});  // 10
        std::vector<std::vector<int>> loops = cubeLoops();
        loops[1] = {4, 5, 8, 9, 10};
        expect(ninja::cellGeometryIsValid(loops, pts), "convex 5-gon boundary face is accepted");
    }

    // 4. An inverted face (wound the wrong way) must fail -- the pyramid test,
    //    kept working alongside the new fan test.
    {
        std::vector<Vec3> pts = cubePoints();
        std::vector<std::vector<int>> loops = cubeLoops();
        std::reverse(loops[1].begin(), loops[1].end());
        expect(!ninja::cellGeometryIsValid(loops, pts), "inward-wound face is REJECTED");
    }

    if (failures == 0) {
        std::printf("cell_geometry: all checks passed\n");
        return 0;
    }
    std::printf("cell_geometry: %d check(s) FAILED\n", failures);
    return 1;
}
