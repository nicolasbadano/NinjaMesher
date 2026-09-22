// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#include "Geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include "CutData.hpp" // Vec3 operators/dot (shared, see CutData.hpp)

namespace ninja {

double pointTriDistSq(const Vec3& p, const Triangle& tri, Vec3* closestOut) {
    const Vec3 ab = tri.v1 - tri.v0;
    const Vec3 ac = tri.v2 - tri.v0;
    const Vec3 ap = p - tri.v0;
    const double d1 = dot(ab, ap);
    const double d2 = dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) {
        if (closestOut) *closestOut = tri.v0;
        const Vec3 d = p - tri.v0;
        return dot(d, d);
    }
    const Vec3 bp = p - tri.v1;
    const double d3 = dot(ab, bp);
    const double d4 = dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) {
        if (closestOut) *closestOut = tri.v1;
        const Vec3 d = p - tri.v1;
        return dot(d, d);
    }
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double v = d1 / (d1 - d3);
        const Vec3 proj = tri.v0 + ab * v;
        if (closestOut) *closestOut = proj;
        const Vec3 d = p - proj;
        return dot(d, d);
    }
    const Vec3 cp = p - tri.v2;
    const double d5 = dot(ab, cp);
    const double d6 = dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) {
        if (closestOut) *closestOut = tri.v2;
        const Vec3 d = p - tri.v2;
        return dot(d, d);
    }
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double w = d2 / (d2 - d6);
        const Vec3 proj = tri.v0 + ac * w;
        if (closestOut) *closestOut = proj;
        const Vec3 d = p - proj;
        return dot(d, d);
    }
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        const Vec3 proj = tri.v1 + (tri.v2 - tri.v1) * w;
        if (closestOut) *closestOut = proj;
        const Vec3 d = p - proj;
        return dot(d, d);
    }
    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom;
    const double w = vc * denom;
    const Vec3 proj = tri.v0 + ab * v + ac * w;
    if (closestOut) *closestOut = proj;
    const Vec3 d = p - proj;
    return dot(d, d);
}

namespace {

// Shared bin-index helper (clamped, same convention throughout this file).
inline int aabbBinIndex1D(double c, double loC, double cellDim, int n) {
    return std::clamp(static_cast<int>((c - loC) / cellDim), 0, n - 1);
}

} // namespace

TriangleAabbBins buildTriangleAabbBins(const std::vector<Triangle>& tris, int targetTrisPerBin) {
    TriangleAabbBins tb;
    tb.tris = &tris;
    if (tris.empty()) {
        tb.nx = tb.ny = tb.nz = 1;
        tb.cellSize = Vec3{1, 1, 1};
        tb.bins.assign(1, {});
        return tb;
    }
    Vec3 lo{1e300, 1e300, 1e300};
    Vec3 hi{-1e300, -1e300, -1e300};
    for (const Triangle& t : tris) {
        for (const Vec3& v : {t.v0, t.v1, t.v2}) {
            lo.x = std::min(lo.x, v.x);
            lo.y = std::min(lo.y, v.y);
            lo.z = std::min(lo.z, v.z);
            hi.x = std::max(hi.x, v.x);
            hi.y = std::max(hi.y, v.y);
            hi.z = std::max(hi.z, v.z);
        }
    }
    // Relative+absolute padding style -- guards
    // against boundary triangles/queries landing exactly on the grid's
    // outer edge (floating point) rather than encoding any tolerance from
    // CutData.cpp's math.
    const Vec3 pad{(hi.x - lo.x) * 1e-6 + 1e-9, (hi.y - lo.y) * 1e-6 + 1e-9, (hi.z - lo.z) * 1e-6 + 1e-9};
    lo = lo - pad;
    hi = hi + pad;
    tb.lo = lo;
    tb.hi = hi;

    const int target = std::max(
        1, static_cast<int>(std::cbrt(static_cast<double>(tris.size()) / std::max(1, targetTrisPerBin))));
    tb.nx = tb.ny = tb.nz = std::max(1, target);
    tb.cellSize = Vec3{(hi.x - lo.x) / tb.nx, (hi.y - lo.y) / tb.ny, (hi.z - lo.z) / tb.nz};
    tb.bins.assign(static_cast<std::size_t>(tb.nx) * static_cast<std::size_t>(tb.ny) * static_cast<std::size_t>(tb.nz),
                    {});

    for (std::size_t i = 0; i < tris.size(); ++i) {
        const Triangle& t = tris[i];
        Vec3 tlo{std::min({t.v0.x, t.v1.x, t.v2.x}), std::min({t.v0.y, t.v1.y, t.v2.y}),
                 std::min({t.v0.z, t.v1.z, t.v2.z})};
        Vec3 thi{std::max({t.v0.x, t.v1.x, t.v2.x}), std::max({t.v0.y, t.v1.y, t.v2.y}),
                 std::max({t.v0.z, t.v1.z, t.v2.z})};
        // Small pad on the triangle's own AABB: guards against a triangle
        // exactly grazing a bin boundary being registered in only one of
        // the two bins that geometrically touch it. Purely a numerical
        // safety margin, same spirit as the grid-level pad above -- it
        // only ever ADDS candidate bins, never removes them, so it cannot
        // break the superset argument.
        const double epx = tb.cellSize.x * 1e-6 + 1e-12;
        const double epy = tb.cellSize.y * 1e-6 + 1e-12;
        const double epz = tb.cellSize.z * 1e-6 + 1e-12;
        const int bx0 = aabbBinIndex1D(tlo.x - epx, lo.x, tb.cellSize.x, tb.nx);
        const int bx1 = aabbBinIndex1D(thi.x + epx, lo.x, tb.cellSize.x, tb.nx);
        const int by0 = aabbBinIndex1D(tlo.y - epy, lo.y, tb.cellSize.y, tb.ny);
        const int by1 = aabbBinIndex1D(thi.y + epy, lo.y, tb.cellSize.y, tb.ny);
        const int bz0 = aabbBinIndex1D(tlo.z - epz, lo.z, tb.cellSize.z, tb.nz);
        const int bz1 = aabbBinIndex1D(thi.z + epz, lo.z, tb.cellSize.z, tb.nz);
        for (int bz = bz0; bz <= bz1; ++bz) {
            for (int by = by0; by <= by1; ++by) {
                for (int bx = bx0; bx <= bx1; ++bx) {
                    const std::size_t idx = static_cast<std::size_t>(bx) + static_cast<std::size_t>(by) * tb.nx +
                                             static_cast<std::size_t>(bz) * tb.nx * tb.ny;
                    tb.bins[idx].push_back(static_cast<int>(i));
                }
            }
        }
    }
    return tb;
}

std::vector<int> queryAabbTriangles(const TriangleAabbBins& bins, const Vec3& loIn, const Vec3& hiIn) {
    std::vector<int> result;
    if (bins.tris == nullptr || bins.tris->empty()) {
        return result;
    }
    // Small pad matching the grid's own numerical margin (see
    // buildTriangleAabbBins) so a query box edge that lands exactly on a
    // bin boundary still visits both adjacent bins.
    const double epx = bins.cellSize.x * 1e-6 + 1e-12;
    const double epy = bins.cellSize.y * 1e-6 + 1e-12;
    const double epz = bins.cellSize.z * 1e-6 + 1e-12;
    const int bx0 = aabbBinIndex1D(loIn.x - epx, bins.lo.x, bins.cellSize.x, bins.nx);
    const int bx1 = aabbBinIndex1D(hiIn.x + epx, bins.lo.x, bins.cellSize.x, bins.nx);
    const int by0 = aabbBinIndex1D(loIn.y - epy, bins.lo.y, bins.cellSize.y, bins.ny);
    const int by1 = aabbBinIndex1D(hiIn.y + epy, bins.lo.y, bins.cellSize.y, bins.ny);
    const int bz0 = aabbBinIndex1D(loIn.z - epz, bins.lo.z, bins.cellSize.z, bins.nz);
    const int bz1 = aabbBinIndex1D(hiIn.z + epz, bins.lo.z, bins.cellSize.z, bins.nz);
    for (int bz = bz0; bz <= bz1; ++bz) {
        for (int by = by0; by <= by1; ++by) {
            for (int bx = bx0; bx <= bx1; ++bx) {
                const std::size_t idx = static_cast<std::size_t>(bx) + static_cast<std::size_t>(by) * bins.nx +
                                         static_cast<std::size_t>(bz) * bins.nx * bins.ny;
                for (int t : bins.bins[idx]) {
                    result.push_back(t);
                }
            }
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

std::vector<int> queryRayTriangles(const TriangleAabbBins& bins, const Vec3& origin, const Vec3& dir) {
    std::vector<int> result;
    if (bins.tris == nullptr || bins.tris->empty()) {
        return result;
    }
    // Slab test: clip the semi-infinite ray (t >= 0) to the bins' overall
    // box. No triangle lies outside this box (every triangle's AABB is
    // inside it by construction), so nothing beyond it can ever be hit.
    double tmin = 0.0;
    double tmax = std::numeric_limits<double>::max();
    const double* o = &origin.x;
    const double* d = &dir.x;
    const double* loP = &bins.lo.x;
    const double* hiP = &bins.hi.x;
    for (int ax = 0; ax < 3; ++ax) {
        if (std::fabs(d[ax]) < 1e-300) {
            if (o[ax] < loP[ax] || o[ax] > hiP[ax]) {
                return result; // parallel to this axis, outside the slab
            }
            continue;
        }
        double t0 = (loP[ax] - o[ax]) / d[ax];
        double t1 = (hiP[ax] - o[ax]) / d[ax];
        if (t0 > t1) std::swap(t0, t1);
        tmin = std::max(tmin, t0);
        tmax = std::min(tmax, t1);
        if (tmin > tmax) {
            return result;
        }
    }
    if (tmax < 0.0) {
        return result;
    }
    tmin = std::max(tmin, 0.0);

    // Entry point, nudged a hair inside the box to avoid landing exactly
    // on a boundary and mis-resolving the starting cell.
    const Vec3 entry{origin.x + dir.x * tmin, origin.y + dir.y * tmin, origin.z + dir.z * tmin};
    int ix = aabbBinIndex1D(entry.x, bins.lo.x, bins.cellSize.x, bins.nx);
    int iy = aabbBinIndex1D(entry.y, bins.lo.y, bins.cellSize.y, bins.ny);
    int iz = aabbBinIndex1D(entry.z, bins.lo.z, bins.cellSize.z, bins.nz);

    auto axisStep = [](double dc, double originC, double loC, double cellDim, int idx, int& step, double& tMaxOut,
                        double& tDeltaOut) {
        if (std::fabs(dc) < 1e-300) {
            step = 0;
            tMaxOut = std::numeric_limits<double>::max();
            tDeltaOut = std::numeric_limits<double>::max();
            return;
        }
        step = dc > 0.0 ? 1 : -1;
        const double boundary = loC + (idx + (dc > 0.0 ? 1 : 0)) * cellDim;
        tMaxOut = (boundary - originC) / dc;
        tDeltaOut = cellDim / std::fabs(dc);
    };
    int stepX, stepY, stepZ;
    double tMaxX, tMaxY, tMaxZ, tDeltaX, tDeltaY, tDeltaZ;
    axisStep(dir.x, origin.x, bins.lo.x, bins.cellSize.x, ix, stepX, tMaxX, tDeltaX);
    axisStep(dir.y, origin.y, bins.lo.y, bins.cellSize.y, iy, stepY, tMaxY, tDeltaY);
    axisStep(dir.z, origin.z, bins.lo.z, bins.cellSize.z, iz, stepZ, tMaxZ, tDeltaZ);

    // Safety cap on step count: the ray can cross at most nx+ny+nz bin
    // boundaries inside the box, plus a small margin for float slop.
    const long maxSteps = static_cast<long>(bins.nx) + bins.ny + bins.nz + 4;
    for (long step = 0; step <= maxSteps; ++step) {
        if (ix < 0 || ix >= bins.nx || iy < 0 || iy >= bins.ny || iz < 0 || iz >= bins.nz) {
            break;
        }
        const std::size_t idx = static_cast<std::size_t>(ix) + static_cast<std::size_t>(iy) * bins.nx +
                                 static_cast<std::size_t>(iz) * bins.nx * bins.ny;
        for (int t : bins.bins[idx]) {
            result.push_back(t);
        }
        // Stop once we've walked past the box exit (tmax).
        const double nextT = std::min({tMaxX, tMaxY, tMaxZ});
        if (nextT > tmax) {
            break;
        }
        if (tMaxX <= tMaxY && tMaxX <= tMaxZ) {
            ix += stepX;
            tMaxX += tDeltaX;
        } else if (tMaxY <= tMaxZ) {
            iy += stepY;
            tMaxY += tDeltaY;
        } else {
            iz += stepZ;
            tMaxZ += tDeltaZ;
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

bool anyTriangleWithin(const TriangleAabbBins& bins, const Vec3& p, double d) {
    if (bins.tris == nullptr || bins.tris->empty()) {
        return false;
    }
    const double dSq = d * d;
    // Same padding convention as queryAabbTriangles: guard against the
    // query cube's edge landing exactly on a bin boundary.
    const double epx = bins.cellSize.x * 1e-6 + 1e-12;
    const double epy = bins.cellSize.y * 1e-6 + 1e-12;
    const double epz = bins.cellSize.z * 1e-6 + 1e-12;
    const Vec3 qlo{p.x - d, p.y - d, p.z - d};
    const Vec3 qhi{p.x + d, p.y + d, p.z + d};
    const int bx0 = aabbBinIndex1D(qlo.x - epx, bins.lo.x, bins.cellSize.x, bins.nx);
    const int bx1 = aabbBinIndex1D(qhi.x + epx, bins.lo.x, bins.cellSize.x, bins.nx);
    const int by0 = aabbBinIndex1D(qlo.y - epy, bins.lo.y, bins.cellSize.y, bins.ny);
    const int by1 = aabbBinIndex1D(qhi.y + epy, bins.lo.y, bins.cellSize.y, bins.ny);
    const int bz0 = aabbBinIndex1D(qlo.z - epz, bins.lo.z, bins.cellSize.z, bins.nz);
    const int bz1 = aabbBinIndex1D(qhi.z + epz, bins.lo.z, bins.cellSize.z, bins.nz);
    for (int bz = bz0; bz <= bz1; ++bz) {
        for (int by = by0; by <= by1; ++by) {
            for (int bx = bx0; bx <= bx1; ++bx) {
                const std::size_t idx = static_cast<std::size_t>(bx) + static_cast<std::size_t>(by) * bins.nx +
                                         static_cast<std::size_t>(bz) * bins.nx * bins.ny;
                for (int t : bins.bins[idx]) {
                    if (pointTriDistSq(p, (*bins.tris)[static_cast<std::size_t>(t)]) <= dSq) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

const std::vector<SphereDir>& sphere14Directions() {
    static const std::vector<SphereDir> dirs = [] {
        std::vector<SphereDir> d;
        d.reserve(14);
        const double faceW = 1.0 / 15.0;
        const double vertW = 3.0 / 40.0;
        for (double s : {1.0, -1.0}) {
            d.push_back({Vec3{s, 0, 0}, faceW});
            d.push_back({Vec3{0, s, 0}, faceW});
            d.push_back({Vec3{0, 0, s}, faceW});
        }
        const double c = 1.0 / std::sqrt(3.0);
        for (double sx : {1.0, -1.0}) {
            for (double sy : {1.0, -1.0}) {
                for (double sz : {1.0, -1.0}) {
                    d.push_back({Vec3{sx * c, sy * c, sz * c}, vertW});
                }
            }
        }
        return d;
    }();
    return dirs;
}

std::vector<double> smoothedSignedDistanceBatch(const TriangleAabbBins& bins, const std::vector<Triangle>& tris,
                                                  const Vec3& locationInMesh, const std::vector<Vec3>& p,
                                                  const std::vector<double>& r, const SmoothedField& cfg,
                                                  std::vector<Vec3>* gradOut) {
    const std::size_t n = p.size();
    std::vector<double> result(n, 0.0);
    if (gradOut) gradOut->assign(n, Vec3{0, 0, 0});
    if (n == 0) return result;

    const std::vector<SphereDir>& dirs = sphere14Directions();
    // Radial quadrature: a simple 3-point rule on [0, r] --
    // centre, r/2 shell, r shell -- with weights 1/6, 1/2, 1/3 (sums to
    // 1). This is Simpson's rule applied to the radial AVERAGING weight
    // (not to phi itself): treating the ball average as a weighted sum
    // over concentric shells at 0, r/2, r with Simpson coefficients
    // 1, 4, 1 over 6 gives exactly this partition, so the two-shell rule
    // integrates a quadratic radial profile of phi exactly, same as
    // Simpson's rule does for a 1-D integrand.
    //
    // test_smoothing also measures the cheaper SINGLE-shell variant
    // (centre + `r` shell only, no `r/2` shell). Its weights are not
    // derived anywhere, so this is a plain two-point split, centre 1/3 / shell
    // 2/3 (documented here as an unforced choice, not derived from any
    // exactness argument beyond "sums to 1 and is antipodally
    // symmetric per direction" -- test_smoothing's measurement decides
    // whether the simplification is acceptable, not this weight choice).
    const double centreW = cfg.singleShell ? (1.0 / 3.0) : (1.0 / 6.0);
    const double halfShellW = cfg.singleShell ? 0.0 : 0.5;
    const double fullShellW = cfg.singleShell ? (2.0 / 3.0) : (1.0 / 3.0);

    // Build the full stencil-point list for every query in one pass, so
    // the sign classification below is ONE batched `classifyVertices`
    // call over every point this evaluation needs -- every stencil point
    // for a whole front goes into one batch, never a second parity test.
    std::vector<Vec3> stencil;
    stencil.reserve(n * (cfg.singleShell ? 15 : 29));
    std::vector<std::size_t> centreIdx(n);
    std::vector<std::array<std::size_t, 14>> halfIdx(cfg.singleShell ? 0 : n);
    std::vector<std::array<std::size_t, 14>> fullIdx(n);
    for (std::size_t i = 0; i < n; ++i) {
        centreIdx[i] = stencil.size();
        stencil.push_back(p[i]);
        if (!cfg.singleShell) {
            for (int k = 0; k < 14; ++k) {
                halfIdx[i][static_cast<std::size_t>(k)] = stencil.size();
                stencil.push_back(p[i] + dirs[static_cast<std::size_t>(k)].o * (0.5 * r[i]));
            }
        }
        for (int k = 0; k < 14; ++k) {
            fullIdx[i][static_cast<std::size_t>(k)] = stencil.size();
            stencil.push_back(p[i] + dirs[static_cast<std::size_t>(k)].o * r[i]);
        }
    }

    // ONE batched parity classification over the whole stencil (the
    // exact ray test, `CutData.cpp`'s `classifyVertices` -- see
    // Geometry.hpp's comment for why this must not be re-implemented).
    const std::vector<bool> stencilSolid = classifyVertices(stencil, tris, locationInMesh);

    // Exact closest-point-on-soup query per stencil point (unbatched --
    // only the SIGN needs batching; the distance query
    // already has no cheaper exact form). `phi = (solid ? -1 : +1) *
    // sqrt(distSq)`.
    // Pure per-point queries; parallel over the stencil (the batch on a big
    // front is hundreds of thousands of points, 29 per query).
    std::vector<double> phi(stencil.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 512)
#endif
    for (std::ptrdiff_t si = 0; si < static_cast<std::ptrdiff_t>(stencil.size()); ++si) {
        const std::size_t i = static_cast<std::size_t>(si);
        const ClosestHit hit = closestPointOnSoup(bins, stencil[i]);
        const double d = std::sqrt(hit.distSq);
        phi[i] = stencilSolid[i] ? -d : d;
    }

    for (std::size_t i = 0; i < n; ++i) {
        double acc = centreW * phi[centreIdx[i]];
        if (!cfg.singleShell) {
            for (int k = 0; k < 14; ++k) {
                acc += halfShellW * dirs[static_cast<std::size_t>(k)].w * phi[halfIdx[i][static_cast<std::size_t>(k)]];
            }
        }
        for (int k = 0; k < 14; ++k) {
            acc += fullShellW * dirs[static_cast<std::size_t>(k)].w * phi[fullIdx[i][static_cast<std::size_t>(k)]];
        }
        result[i] = acc;

        if (gradOut && r[i] > 1e-300) {
            // S = sum_k w_k * o_k * phi(p + r*o_k) over the FULL (`r`)
            // shell only -- see the header comment for the 3*S/r
            // derivation.
            Vec3 s{0, 0, 0};
            for (int k = 0; k < 14; ++k) {
                const SphereDir& d = dirs[static_cast<std::size_t>(k)];
                s = s + d.o * (d.w * phi[fullIdx[i][static_cast<std::size_t>(k)]]);
            }
            (*gradOut)[i] = s * (3.0 / r[i]);
        }
    }
    return result;
}

double smoothedSignedDistance(const TriangleAabbBins& bins, const std::vector<Triangle>& tris,
                                const Vec3& locationInMesh, const Vec3& p, double r, const SmoothedField& cfg) {
    const std::vector<Vec3> pv{p};
    const std::vector<double> rv{r};
    return smoothedSignedDistanceBatch(bins, tris, locationInMesh, pv, rv, cfg)[0];
}

ClosestHit closestPointOnSoup(const TriangleAabbBins& bins, const Vec3& p) {
    ClosestHit best;
    if (bins.tris == nullptr || bins.tris->empty()) {
        return best;
    }

    const int bx = aabbBinIndex1D(p.x, bins.lo.x, bins.cellSize.x, bins.nx);
    const int by = aabbBinIndex1D(p.y, bins.lo.y, bins.cellSize.y, bins.ny);
    const int bz = aabbBinIndex1D(p.z, bins.lo.z, bins.cellSize.z, bins.nz);

    auto testBin = [&](int ix, int iy, int iz) {
        if (ix < 0 || ix >= bins.nx || iy < 0 || iy >= bins.ny || iz < 0 || iz >= bins.nz) {
            return;
        }
        const std::size_t idx = static_cast<std::size_t>(ix) + static_cast<std::size_t>(iy) * bins.nx +
                                 static_cast<std::size_t>(iz) * bins.nx * bins.ny;
        for (int t : bins.bins[idx]) {
            Vec3 closest;
            const double dSq = pointTriDistSq(p, (*bins.tris)[static_cast<std::size_t>(t)], &closest);
            if (dSq < best.distSq || (dSq == best.distSq && t < best.triangle)) {
                best.distSq = dSq;
                best.point = closest;
                best.triangle = t;
            }
        }
    };

    // Ring 0: the query point's own (clamped) bin.
    testBin(bx, by, bz);

    int prevX0 = bx, prevX1 = bx, prevY0 = by, prevY1 = by, prevZ0 = bz, prevZ1 = bz;
    const int maxRadius = std::max({bins.nx, bins.ny, bins.nz});
    for (int r = 1; r <= maxRadius; ++r) {
        const int x0 = std::max(0, bx - r), x1 = std::min(bins.nx - 1, bx + r);
        const int y0 = std::max(0, by - r), y1 = std::min(bins.ny - 1, by + r);
        const int z0 = std::max(0, bz - r), z1 = std::min(bins.nz - 1, bz + r);

        // Visit only the NEW bins in this ring (box(r) minus box(r-1)) --
        // iterate the full box and skip anything already inside the
        // previous (smaller) box, which is cheap at this project's bin
        // counts and keeps the logic simple/obviously-correct.
        for (int iz = z0; iz <= z1; ++iz) {
            for (int iy = y0; iy <= y1; ++iy) {
                for (int ix = x0; ix <= x1; ++ix) {
                    const bool wasVisited = ix >= prevX0 && ix <= prevX1 && iy >= prevY0 && iy <= prevY1 &&
                                             iz >= prevZ0 && iz <= prevZ1;
                    if (wasVisited) {
                        continue;
                    }
                    testBin(ix, iy, iz);
                }
            }
        }
        prevX0 = x0; prevX1 = x1; prevY0 = y0; prevY1 = y1; prevZ0 = z0; prevZ1 = z1;

        // Pruning check (see the correctness argument in Geometry.hpp):
        // the world-space box covered by bins visited so far.
        const Vec3 boxLo{bins.lo.x + x0 * bins.cellSize.x, bins.lo.y + y0 * bins.cellSize.y,
                          bins.lo.z + z0 * bins.cellSize.z};
        const Vec3 boxHi{bins.lo.x + (x1 + 1) * bins.cellSize.x, bins.lo.y + (y1 + 1) * bins.cellSize.y,
                          bins.lo.z + (z1 + 1) * bins.cellSize.z};
        const double boundaryDist = std::min({p.x - boxLo.x, boxHi.x - p.x, p.y - boxLo.y, boxHi.y - p.y,
                                               p.z - boxLo.z, boxHi.z - p.z});
        const bool boxCoversGrid = x0 == 0 && x1 == bins.nx - 1 && y0 == 0 && y1 == bins.ny - 1 && z0 == 0 &&
                                    z1 == bins.nz - 1;
        if (boxCoversGrid) {
            break; // nothing left to visit
        }
        if (best.triangle >= 0 && boundaryDist > 0.0 && boundaryDist * boundaryDist >= best.distSq) {
            break; // no unvisited bin can hold a closer triangle
        }
    }
    return best;
}

} // namespace ninja
