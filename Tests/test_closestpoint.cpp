// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

// Unit-style ctest for closestPointOnSoup.
//
// Usage: test_closestpoint <icosphere.stl path>
//
// Part A: hand-computable point/triangle configurations, asserted
// exactly (face-interior, edge, vertex, equidistant tie).
// Part B: brute-force comparison over the icosphere fixture at 1000
// deterministic-seeded random points -- bitwise-equal distSq and
// triangle index vs. a plain O(n) scan with the SAME lowest-index
// tie-break rule.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>

#include "CutData.hpp" // Vec3 operators
#include "Geometry.hpp"
#include "Stl.hpp"

using namespace ninja;

namespace {

int failures = 0;

void expectNear(double a, double b, double tol, const char* what) {
    if (std::fabs(a - b) > tol) {
        std::fprintf(stderr, "FAIL: %s: got %.17g expected %.17g (tol %.3g)\n", what, a, b, tol);
        ++failures;
    }
}

void expectEq(int a, int b, const char* what) {
    if (a != b) {
        std::fprintf(stderr, "FAIL: %s: got %d expected %d\n", what, a, b);
        ++failures;
    }
}

// Brute-force reference: lowest-triangle-index tie-break, identical rule
// to closestPointOnSoup's.
ClosestHit bruteForce(const std::vector<Triangle>& tris, const Vec3& p) {
    ClosestHit best;
    for (std::size_t i = 0; i < tris.size(); ++i) {
        Vec3 closest;
        const double dSq = pointTriDistSq(p, tris[i], &closest);
        if (dSq < best.distSq || (dSq == best.distSq && static_cast<int>(i) < best.triangle)) {
            best.distSq = dSq;
            best.point = closest;
            best.triangle = static_cast<int>(i);
        }
    }
    return best;
}

void partA() {
    // A single unit triangle in the z=0 plane: (0,0,0)-(1,0,0)-(0,1,0).
    std::vector<Triangle> tris(1);
    tris[0].v0 = Vec3{0, 0, 0};
    tris[0].v1 = Vec3{1, 0, 0};
    tris[0].v2 = Vec3{0, 1, 0};
    tris[0].solidId = 0;
    TriangleAabbBins bins = buildTriangleAabbBins(tris, 1);

    // Face-interior: point above the centroid.
    {
        ClosestHit h = closestPointOnSoup(bins, Vec3{0.2, 0.2, 1.0});
        expectNear(h.distSq, 1.0, 1e-12, "face-interior distSq");
        expectNear(h.point.x, 0.2, 1e-12, "face-interior point.x");
        expectNear(h.point.y, 0.2, 1e-12, "face-interior point.y");
        expectNear(h.point.z, 0.0, 1e-12, "face-interior point.z");
        expectEq(h.triangle, 0, "face-interior triangle");
    }
    // Edge: point above the midpoint of the hypotenuse (v1-v2), offset
    // beyond the triangle along the outward edge normal.
    {
        // hypotenuse midpoint (0.5,0.5,0), outward in-plane normal
        // direction is (1,1,0)/sqrt(2).
        const double s = 1.0 / std::sqrt(2.0);
        Vec3 q{0.5 + 0.3 * s, 0.5 + 0.3 * s, 0.0};
        ClosestHit h = closestPointOnSoup(bins, q);
        expectNear(h.distSq, 0.09, 1e-9, "edge distSq");
        expectNear(h.point.x, 0.5, 1e-9, "edge point.x");
        expectNear(h.point.y, 0.5, 1e-9, "edge point.y");
        expectEq(h.triangle, 0, "edge triangle");
    }
    // Vertex: point beyond v1=(1,0,0), along (1,0,0) extended.
    {
        Vec3 q{1.5, -0.5, 0.0}; // closest feature is the vertex v1
        ClosestHit h = closestPointOnSoup(bins, q);
        expectNear(h.point.x, 1.0, 1e-12, "vertex point.x");
        expectNear(h.point.y, 0.0, 1e-12, "vertex point.y");
        expectNear(h.distSq, 0.25 + 0.25, 1e-12, "vertex distSq");
        expectEq(h.triangle, 0, "vertex triangle");
    }
    // Equidistant tie: two triangles, symmetric about the query point,
    // exactly equal distance -- lowest index must win. Triangle 0 at
    // x=-2 (a small triangle whose closest point is its vertex at
    // x=-1), triangle 1 mirrored at x=+1.
    {
        std::vector<Triangle> tt(2);
        tt[0].v0 = Vec3{-1, 0, 0};
        tt[0].v1 = Vec3{-2, 1, 0};
        tt[0].v2 = Vec3{-2, -1, 0};
        tt[1].v0 = Vec3{1, 0, 0};
        tt[1].v1 = Vec3{2, 1, 0};
        tt[1].v2 = Vec3{2, -1, 0};
        TriangleAabbBins tb = buildTriangleAabbBins(tt, 1);
        ClosestHit h = closestPointOnSoup(tb, Vec3{0, 0, 0});
        expectEq(h.triangle, 0, "tie-break lowest index");
        expectNear(h.distSq, 1.0, 1e-12, "tie distSq");

        // And the reverse declaration order: swap which one is index 0,
        // confirming the tie-break tracks the INDEX, not geometry.
        std::vector<Triangle> tt2(2);
        tt2[0] = tt[1];
        tt2[1] = tt[0];
        TriangleAabbBins tb2 = buildTriangleAabbBins(tt2, 1);
        ClosestHit h2 = closestPointOnSoup(tb2, Vec3{0, 0, 0});
        expectEq(h2.triangle, 0, "tie-break lowest index (swapped)");
    }
}

void partB(const std::string& stlPath) {
    std::vector<Triangle> tris = readStl(stlPath);
    if (tris.empty()) {
        std::fprintf(stderr, "FAIL: icosphere fixture '%s' has no triangles\n", stlPath.c_str());
        ++failures;
        return;
    }
    TriangleAabbBins bins = buildTriangleAabbBins(tris);

    std::mt19937 rng(1234567u); // fixed seed -- deterministic
    std::uniform_real_distribution<double> coord(-0.6, 1.6); // spans inside/near/outside r=0.3 sphere at (0.5,.5,.5)

    constexpr int kN = 1000;
    for (int i = 0; i < kN; ++i) {
        Vec3 p{coord(rng), coord(rng), coord(rng)};
        ClosestHit fast = closestPointOnSoup(bins, p);
        ClosestHit ref = bruteForce(tris, p);
        if (fast.distSq != ref.distSq || fast.triangle != ref.triangle) {
            std::fprintf(stderr,
                         "FAIL: sample %d p=(%.17g %.17g %.17g) fast(distSq=%.17g tri=%d) != "
                         "brute(distSq=%.17g tri=%d)\n",
                         i, p.x, p.y, p.z, fast.distSq, fast.triangle, ref.distSq, ref.triangle);
            ++failures;
        }
        if (fast.point.x != ref.point.x || fast.point.y != ref.point.y || fast.point.z != ref.point.z) {
            std::fprintf(stderr, "FAIL: sample %d point mismatch fast=(%.17g %.17g %.17g) ref=(%.17g %.17g %.17g)\n",
                         i, fast.point.x, fast.point.y, fast.point.z, ref.point.x, ref.point.y, ref.point.z);
            ++failures;
        }
    }
    if (failures == 0) {
        std::printf("PASS: %d brute-force-matched random samples\n", kN);
    }
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: test_closestpoint <icosphere.stl>\n");
        return 1;
    }
    partA();
    partB(argv[1]);
    if (failures == 0) {
        std::printf("PASS: all hand-computed configurations exact\n");
        return 0;
    }
    std::fprintf(stderr, "%d failure(s)\n", failures);
    return 1;
}
