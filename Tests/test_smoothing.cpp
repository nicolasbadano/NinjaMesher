// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

// The smoothed signed-distance field, STANDALONE and
// MEASURED before any wiring into the march (`src/Layers.cpp`) --
// mirrors `test_closestpoint.cpp` for the closest-point query
// and the design's own "measure-first" gate (the design predates the
// smoothing rename, so its wording differs).
//
// Usage: test_smoothing <icosphere.stl path>   (sphere radius 0.3 @
//                                              (0.5,0.5,0.5), 3 subdiv --
//                                              same fixture as
//                                              test_closestpoint)
//
// Cases:
//  A. flat wall: phi_r == phi to round-off (the single most important
//     test -- a smoothed flat wall must not move).
//  B. sphere of radius R: phi_r recovers phi to O(r^2/R); the measured
//     constant is printed, not just checked against a bound.
//  C. 90-degree convex edge: corner displacement ~ r*(sqrt2-1).
//  D. tessellation invariance: 2-triangle vs ~2000-triangle flat wall,
//     same phi_r to round-off (the test that would have caught the
//     soft-min trap the design rejects).
//  E. sign correctness straddling the surface.
// Plus the Task-0-mandated report: per-evaluation cost relative to a bare
// closestPointOnSoup call, and the single-shell-vs-two-shell comparison
// on cases B/C above.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <random>
#include <vector>

#include "CutData.hpp" // Vec3 operators
#include "Geometry.hpp"
#include "Stl.hpp"

using namespace ninja;

namespace {

int failures = 0;

void expectNear(double a, double b, double tol, const char* what) {
    if (std::fabs(a - b) > tol) {
        std::fprintf(stderr, "FAIL: %s: got %.17g expected %.17g (tol %.3g, diff %.3g)\n", what, a, b, tol,
                     std::fabs(a - b));
        ++failures;
    }
}

void expectTrue(bool cond, const char* what) {
    if (!cond) {
        std::fprintf(stderr, "FAIL: %s\n", what);
        ++failures;
    }
}

// Axis-aligned box, 12 triangles (2/face), same construction as
// Tests/scripts/make_stls.py's box_triangles -- a closed, watertight
// solid so classifyVertices' ray-parity test is well-defined.
std::vector<Triangle> boxTriangles(Vec3 lo, Vec3 hi) {
    Vec3 p000{lo.x, lo.y, lo.z}, p100{hi.x, lo.y, lo.z}, p010{lo.x, hi.y, lo.z}, p110{hi.x, hi.y, lo.z};
    Vec3 p001{lo.x, lo.y, hi.z}, p101{hi.x, lo.y, hi.z}, p011{lo.x, hi.y, hi.z}, p111{hi.x, hi.y, hi.z};
    auto quad = [](std::vector<Triangle>& out, Vec3 a, Vec3 b, Vec3 c, Vec3 d) {
        Triangle t1;
        t1.v0 = a;
        t1.v1 = b;
        t1.v2 = c;
        out.push_back(t1);
        Triangle t2;
        t2.v0 = a;
        t2.v1 = c;
        t2.v2 = d;
        out.push_back(t2);
    };
    std::vector<Triangle> tris;
    quad(tris, p000, p010, p011, p001); // -x
    quad(tris, p100, p101, p111, p110); // +x
    quad(tris, p000, p001, p101, p100); // -y
    quad(tris, p010, p110, p111, p011); // +y
    quad(tris, p000, p100, p110, p010); // -z
    quad(tris, p001, p011, p111, p101); // +z
    for (std::size_t i = 0; i < tris.size(); ++i) tris[i].solidId = 0;
    return tris;
}

// Same +z face as `boxTriangles`, but that ONE face subdivided into an
// nxn grid (~2*n*n triangles) instead of 2 -- everything else identical.
// Tests that phi_r on a flat wall is INDEPENDENT of how finely the wall
// happens to be tessellated (the tessellation-invariance
// requirement; the property a soft-min-over-triangles formulation would
// have failed).
std::vector<Triangle> boxTrianglesFinePlusZ(Vec3 lo, Vec3 hi, int n) {
    std::vector<Triangle> tris = boxTriangles(lo, hi);
    // Drop the 2 coarse +z triangles (indices 10, 11 in the ordering
    // above: quads pushed in order -x,+x,-y,+y,-z,+z, 2 tris each).
    tris.resize(10);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            const double x0 = lo.x + (hi.x - lo.x) * i / n;
            const double x1 = lo.x + (hi.x - lo.x) * (i + 1) / n;
            const double y0 = lo.y + (hi.y - lo.y) * j / n;
            const double y1 = lo.y + (hi.y - lo.y) * (j + 1) / n;
            Vec3 a{x0, y0, hi.z}, b{x1, y0, hi.z}, c{x1, y1, hi.z}, d{x0, y1, hi.z};
            Triangle t1;
            t1.v0 = a;
            t1.v1 = b;
            t1.v2 = c;
            t1.solidId = 0;
            Triangle t2;
            t2.v0 = a;
            t2.v1 = c;
            t2.v2 = d;
            t2.solidId = 0;
            tris.push_back(t1);
            tris.push_back(t2);
        }
    }
    return tris;
}

void caseA_flatWall() {
    // A cube [0,10]^3, evaluated near the centre of the +z face
    // (z=10, x=y=5) -- 5 units from every edge, far beyond any `r`
    // tested here, so the field near this point is exactly the plane
    // z=10 with no corner in range of the stencil.
    const std::vector<Triangle> tris = boxTriangles(Vec3{0, 0, 0}, Vec3{10, 10, 10});
    const TriangleAabbBins bins = buildTriangleAabbBins(tris);
    const Vec3 loc{-1, -1, -1}; // known-fluid reference point, OUTSIDE the box
    const SmoothedField cfgTwo{0.25, false};
    const SmoothedField cfgOne{0.25, true};
    // The box occupies z in [0, 10] (SOLID); fluid is z > 10, so at
    // p = (5, 5, 10 + off) the EXACT signed distance is exactly `off`
    // (>0 fluid, <0 solid) near the flat +z face, away from any edge.
    for (double off : {0.5, 0.1, -0.5, -0.1}) {
        const Vec3 p{5, 5, 10 + off};
        const double r = 0.2; // arbitrary, well inside the 5-unit edge clearance
        const double phiTwo = smoothedSignedDistance(bins, tris, loc, p, r, cfgTwo);
        const double phiOne = smoothedSignedDistance(bins, tris, loc, p, r, cfgOne);
        expectNear(phiTwo, off, 1e-9, "flat wall two-shell phi_r == phi");
        expectNear(phiOne, off, 1e-9, "flat wall single-shell phi_r == phi");
    }
    std::printf("PASS: case A (flat wall, phi_r == phi to round-off, both shell variants)\n");
}

void caseB_sphere(const std::string& sphereStl) {
    const std::vector<Triangle> tris = readStl(sphereStl);
    if (tris.empty()) {
        std::fprintf(stderr, "FAIL: sphere fixture '%s' has no triangles\n", sphereStl.c_str());
        ++failures;
        return;
    }
    const TriangleAabbBins bins = buildTriangleAabbBins(tris);
    const Vec3 centre{0.5, 0.5, 0.5};
    const Vec3 loc{0, 0, 0}; // known-fluid reference, OUTSIDE the r=0.3 sphere
    const double R = 0.3;
    const SmoothedField cfgTwo{0.25, false};
    const SmoothedField cfgOne{0.25, true};
    // Sample several points along +x outside the sphere, at varying `r`,
    // and fit the measured constant in phi_r - phi ~= k * r^2 / R.
    for (double r : {0.01, 0.02, 0.04}) {
        const Vec3 p{centre.x + R + 0.05, centre.y, centre.z}; // exact phi = 0.05
        const double phiExact = 0.05;
        const double phiTwo = smoothedSignedDistance(bins, tris, loc, p, r, cfgTwo);
        const double phiOne = smoothedSignedDistance(bins, tris, loc, p, r, cfgOne);
        const double kTwo = (phiTwo - phiExact) * R / (r * r);
        const double kOne = (phiOne - phiExact) * R / (r * r);
        std::printf(
            "MEASURED case B (sphere R=%.2g): r=%.3g two-shell phi_r=%.6g (k=%.4g)  single-shell "
            "phi_r=%.6g (k=%.4g)  exact=%.6g\n",
            R, r, phiTwo, kTwo, phiOne, kOne, phiExact);
        // The convex sphere is approached from OUTSIDE by a stencil that
        // dips partly inside the sphere at large enough r/R -- smoothing
        // a convex bulge always raises phi_r above the exact value there
        // (the ball average sees more "further away" solid than fluid on
        // net... conversely near a convex bump approached from outside,
        // the far side of the stencil samples the true surface's
        // curvature and phi_r < phi is also plausible depending on
        // sign convention -- reported empirically above rather than
        // asserted a priori beyond a loose sanity bound).
        expectTrue(std::fabs(phiTwo - phiExact) < 0.5 * r, "sphere case: two-shell displacement bounded by 0.5r");
    }
}

void caseC_convexEdge() {
    // Cube corner at the origin, evaluated ALONG the (1,1,1)/sqrt3
    // bisector of the corner (0,0,0)-(0,0,0)... i.e. the corner of
    // [-1,1]^3-ish. Use box [0,2]^3 so the corner of interest is
    // (0,0,0) and its three faces meet there; evaluate along the
    // outward corner bisector direction (-1,-1,-1)/sqrt3 from the
    // corner, in the fluid.
    const std::vector<Triangle> tris = boxTriangles(Vec3{0, 0, 0}, Vec3{2, 2, 2});
    const TriangleAabbBins bins = buildTriangleAabbBins(tris);
    const Vec3 loc{-1, -1, -1}; // known-fluid reference, OUTSIDE the box
    const double r = 0.05;
    const SmoothedField cfgTwo{0.25, false};
    const double c = 1.0 / std::sqrt(3.0);
    // Point sitting exactly on the true corner (the exact d=0 level
    // set's apex): phi there is 0 exactly (the corner point itself).
    const Vec3 corner{0, 0, 0};
    const double phiAtCorner = smoothedSignedDistance(bins, tris, loc, corner, r, cfgTwo);
    // Rounding chamfers the corner AWAY from the solid: the true corner
    // point (exact phi = 0, sitting exactly on the surface) now reads as
    // FLUID (phi_r > 0) by roughly r*(sqrt2-1) (2-plane edge) to
    // r*(sqrt3-1) (3-plane corner) -- the solid shrinks at convex
    // features, exactly the
    // "convex edges round with radius ~r" property the design predicts.
    std::printf("MEASURED case C (90deg/3-plane corner): phi_r(corner) = %.6g at r=%.4g (expected ~ +%.4g to +%.4g)\n",
                 phiAtCorner, r, r * (std::sqrt(2.0) - 1.0), r * (std::sqrt(3.0) - 1.0));
    expectTrue(phiAtCorner > 0.0, "corner: smoothed field reads corner as fluid (rounding), sign positive");
    expectTrue(std::fabs(phiAtCorner) < 1.5 * r * (std::sqrt(3.0) - 1.0),
               "corner displacement roughly matches the r*(sqrt3-1) bound (with slack)");
    (void)c;
}

void caseD_tessellationInvariance() {
    const std::vector<Triangle> coarse = boxTriangles(Vec3{0, 0, 0}, Vec3{10, 10, 10}); // 2 tris on +z
    const std::vector<Triangle> fine = boxTrianglesFinePlusZ(Vec3{0, 0, 0}, Vec3{10, 10, 10}, 32); // ~2048 tris on +z
    const TriangleAabbBins binsCoarse = buildTriangleAabbBins(coarse);
    const TriangleAabbBins binsFine = buildTriangleAabbBins(fine);
    const Vec3 loc{-1, -1, -1}; // known-fluid reference, OUTSIDE the box
    const SmoothedField cfgTwo{0.25, false};
    for (double signedDist : {0.3, 0.05, -0.05}) {
        const Vec3 p{5, 5, 10 - signedDist};
        const double r = 0.15;
        const double phiCoarse = smoothedSignedDistance(binsCoarse, coarse, loc, p, r, cfgTwo);
        const double phiFine = smoothedSignedDistance(binsFine, fine, loc, p, r, cfgTwo);
        std::printf("MEASURED case D (tessellation invariance): 2-tri phi_r=%.15g  ~2000-tri phi_r=%.15g  diff=%.3g\n",
                     phiCoarse, phiFine, std::fabs(phiCoarse - phiFine));
        expectNear(phiCoarse, phiFine, 1e-9, "tessellation invariance: coarse vs fine +z face phi_r");
    }
}

void caseE_signStraddling() {
    const std::vector<Triangle> tris = boxTriangles(Vec3{0, 0, 0}, Vec3{10, 10, 10});
    const TriangleAabbBins bins = buildTriangleAabbBins(tris);
    const Vec3 loc{-1, -1, -1}; // known-fluid reference, OUTSIDE the box
    const SmoothedField cfgTwo{0.25, false};
    const double fluid = smoothedSignedDistance(bins, tris, loc, Vec3{5, 5, 10.2}, 0.05, cfgTwo);
    const double solid = smoothedSignedDistance(bins, tris, loc, Vec3{5, 5, 9.8}, 0.05, cfgTwo);
    expectTrue(fluid > 0.0, "sign: fluid side positive");
    expectTrue(solid < 0.0, "sign: solid side negative");
    std::printf("PASS: case E (sign straddling: fluid=%.6g solid=%.6g)\n", fluid, solid);
}

void costReport(const std::string& sphereStl) {
    const std::vector<Triangle> tris = readStl(sphereStl);
    const TriangleAabbBins bins = buildTriangleAabbBins(tris);
    const Vec3 loc{0, 0, 0}; // known-fluid reference, OUTSIDE the r=0.3 sphere
    constexpr int kN = 2000;
    std::vector<Vec3> pts(kN);
    std::vector<double> radii(kN, 0.02);
    std::mt19937 rng(99u);
    std::uniform_real_distribution<double> coord(0.15, 0.85);
    for (int i = 0; i < kN; ++i) pts[static_cast<std::size_t>(i)] = Vec3{coord(rng), coord(rng), coord(rng)};

    auto t0 = std::chrono::steady_clock::now();
    for (const Vec3& p : pts) { volatile double d = std::sqrt(closestPointOnSoup(bins, p).distSq); (void)d; }
    auto t1 = std::chrono::steady_clock::now();
    const SmoothedField cfgOne{0.25, true};
    const std::vector<double> phiOne = smoothedSignedDistanceBatch(bins, tris, loc, pts, radii, cfgOne);
    auto t2 = std::chrono::steady_clock::now();
    const SmoothedField cfgTwo{0.25, false};
    const std::vector<double> phiTwo = smoothedSignedDistanceBatch(bins, tris, loc, pts, radii, cfgTwo);
    auto t3 = std::chrono::steady_clock::now();

    const double bareUs = std::chrono::duration<double, std::micro>(t1 - t0).count() / kN;
    const double oneUs = std::chrono::duration<double, std::micro>(t2 - t1).count() / kN;
    const double twoUs = std::chrono::duration<double, std::micro>(t3 - t2).count() / kN;
    std::printf(
        "COST (per evaluation, N=%d, icosphere fixture): bare closestPointOnSoup %.3g us; "
        "single-shell (15 pts) %.3g us (%.2fx); two-shell (29 pts) %.3g us (%.2fx)\n",
        kN, bareUs, oneUs, oneUs / bareUs, twoUs, twoUs / bareUs);
    (void)phiOne;
    (void)phiTwo;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: test_smoothing <icosphere.stl>\n");
        return 1;
    }
    // A: quadrature sums to unity (both direction-weight subtotals, plan
    // the quadrature's own stated check).
    {
        const std::vector<SphereDir>& dirs = sphere14Directions();
        double sum = 0.0;
        for (const SphereDir& d : dirs) sum += d.w;
        expectNear(sum, 1.0, 1e-15, "14-point quadrature weights sum to 1");
        expectNear(static_cast<double>(dirs.size()), 14.0, 0.0, "14 directions");
    }
    caseA_flatWall();
    caseB_sphere(argv[1]);
    caseC_convexEdge();
    caseD_tessellationInvariance();
    caseE_signStraddling();
    costReport(argv[1]);

    if (failures == 0) {
        std::printf("PASS: all smoothing checks\n");
        return 0;
    }
    std::fprintf(stderr, "%d failure(s)\n", failures);
    return 1;
}
