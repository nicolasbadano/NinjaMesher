// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <limits>
#include <vector>

#include "BaseMesh.hpp"
#include "Stl.hpp"

namespace ninja {

// Closest squared distance from `p` to triangle `tri` (standard clamped-
// barycentric closest point; robust for degenerate-ish but non-zero-area
// triangles). Shared by Cutter.cpp and the distance-mode refinement
// predicate (Refine.cpp) -- both consumers need "closest point on a soup".
//
// Optional `closestOut`
// output parameter -- when non-null, the closest point on the triangle is
// also written there. Existing callers (Refine.cpp's distance predicate,
// anyTriangleWithin below) pass nullptr and pay no extra cost beyond a
// null check per branch; only closestPointOnSoup's per-triangle test
// below asks for the point.
double pointTriDistSq(const Vec3& p, const Triangle& tri, Vec3* closestOut = nullptr);

// --- AABB-based triangle binning (CutData.cpp acceleration; also used by
// Refine.cpp's distance-mode refinement predicate, see anyTriangleWithin
// below) ---------------------------------------------------------------
//
// A centroid-only bin grid is not exact for point/segment/ray queries: a
// large triangle's centroid can sit in a bin far from a part of its AABB
// that a query actually overlaps, silently dropping a triangle the
// brute-force loop would have tested. That would be tolerable for an
// approximate nearest-distance search, but not for an EXACT predicate, or
// for CutData.cpp's classification/intercept math, which must be
// bit-identical to the old full-soup scan.
//
// TriangleAabbBins fixes this by registering each triangle into EVERY
// bin its full (padded) AABB overlaps, so a query that gathers every bin
// overlapping its own (padded) region is guaranteed to find every
// triangle whose AABB could possibly intersect it — a strict superset of
// the brute-force candidate set, never a subset. The bin grid is a "dumb
// uniform grid" (no kd-tree/BVH).
struct TriangleAabbBins {
    Vec3 lo{0, 0, 0};
    Vec3 hi{0, 0, 0};
    int nx = 1;
    int ny = 1;
    int nz = 1;
    Vec3 cellSize{1, 1, 1};
    std::vector<std::vector<int>> bins; // triangle indices, size nx*ny*nz
    const std::vector<Triangle>* tris = nullptr;
};

// Builds an AABB-based bin grid over `tris`' bounding box, aiming for
// roughly `targetTrisPerBin` triangles per bin (AABB-overlap insertion
// means large triangles land in more than one bin, which is expected and
// harmless -- correctness, not bin-occupancy balance, is the goal).
// `tris` must outlive the returned TriangleAabbBins.
TriangleAabbBins buildTriangleAabbBins(const std::vector<Triangle>& tris, int targetTrisPerBin = 8);

// Returns the deduplicated (sorted, unique) indices of every triangle
// whose (padded) AABB overlaps the query box [lo, hi]. Exact superset
// argument: if a triangle's AABB overlaps [lo, hi], the (non-empty)
// overlap region is covered by the bin grid, so at least one bin
// overlaps both the triangle's AABB (where the triangle was inserted)
// and the query box [lo, hi] (which this function visits) -- the
// triangle is therefore found. Safe to call with lo==hi for a point
// query.
std::vector<int> queryAabbTriangles(const TriangleAabbBins& bins, const Vec3& lo, const Vec3& hi);

// Returns the deduplicated indices of every triangle whose bin(s) are
// crossed by the ray `origin + t*dir`, `t >= 0` (a semi-infinite ray, as
// used by CutData.cpp's parity test). The ray is clipped to the bins'
// overall (padded) bounding box -- no triangle lies outside it, so
// nothing beyond the box can be hit -- then walked bin-to-bin with a
// standard 3D grid (Amanatides-Woo style) traversal from the entry point
// to the exit point, collecting every bin touched. Superset argument:
// any triangle the ray could geometrically hit lies inside the overall
// bounding box, so its true intersection point (if any) lies in some bin
// on the ray's path through the box; that triangle was inserted into
// every bin its AABB overlaps, which includes the bin containing the
// intersection point (or a bin adjacent to it if the hit is exactly on a
// bin boundary -- covered by the same padding used for AABB insertion),
// so it is visited. Returns empty if the ray misses the box entirely.
std::vector<int> queryRayTriangles(const TriangleAabbBins& bins, const Vec3& origin, const Vec3& dir);

// Exact early-exit within-distance predicate (serving Refine.cpp's
// surface-rule marking): true iff at least one triangle in
// `bins` is within `d` of `p` (i.e. exists a triangle with
// pointTriDistSq(p, tri) <= d*d). Only a yes/no answer is needed, so this
// can afford to (a) use the exact AABB-based bins (`TriangleAabbBins`,
// not a centroid-only grid) and (b) stop at the FIRST qualifying triangle
// instead of resolving the true minimum distance.
//
// Exactness argument: any triangle within `d` of `p` has an AABB that
// intersects the axis-aligned cube [p-d, p+d] (the triangle has a point at
// distance <= d from p, and that point lies on/in the triangle's own
// AABB). Visiting every bin whose extent overlaps that cube is therefore a
// superset of "bins containing a triangle within d of p" -- the same
// superset argument as `queryAabbTriangles`. The per-triangle test
// (`pointTriDistSq(p, tri) <= d*d`) resolves the same truth value as a
// full nearest-distance query would -- only less work is done (no ring
// expansion, no running minimum, early exit on the first hit, exact bins
// instead of centroid-only bins).
bool anyTriangleWithin(const TriangleAabbBins& bins, const Vec3& p, double d);

// Exact closest-point-on-soup query. `anyTriangleWithin` only answers a
// boolean; the layers march needs the richer
// {distance, point, owning triangle} answer, hence this separate query.
struct ClosestHit {
    double distSq = std::numeric_limits<double>::max();
    Vec3 point{0, 0, 0};
    int triangle = -1;
};

// Exact nearest point (and owning triangle) on the soup registered in
// `bins`, to query point `p`. Deterministic tie-break: among triangles at
// exactly equal distSq, the LOWEST triangle index wins (checked via
// explicit index comparison in the scan below, independent of bin visit
// order).
//
// Algorithm: expanding-ring search over the AABB bin grid, seeded by the
// bin size. Ring r (r=0,1,2,...) visits exactly the shell of bins whose
// clamped-to-grid index box grew relative to ring r-1 (ring 0 is the
// query point's own bin). After each ring, every triangle visited so far
// is a candidate; a further ring may only be explored if it could still
// contain a strictly closer triangle.
//
// Pruning-correctness argument: after ring r, all visited bins together
// cover a world-space box B (the clamped bin-index box converted to
// coordinates, using the SAME bin grid `buildTriangleAabbBins` used to
// register every triangle by full AABB-overlap, per its own superset
// argument above). Any triangle registered in a bin OUTSIDE B has its
// AABB overlapping ONLY bins outside B (by the registration rule), so no
// point of that triangle can lie inside B's interior at a distance from
// `p` less than the distance from `p` to B's own boundary -- because
// getting from `p` (which is inside B, since B always contains p's own
// bin) to any point outside B requires crossing B's boundary along at
// least one axis, whose distance is exactly
// min(p.x - B.lo.x, B.hi.x - p.x, p.y - B.lo.y, B.hi.y - p.y,
//     p.z - B.lo.z, B.hi.z - p.z).
// If that minimum boundary distance, squared, is >= the current best
// distSq, every unvisited triangle (whose AABB lies entirely outside B)
// is guaranteed no closer than the current best, so the search may stop.
// This is the standard "distance to unexplored region" ring-pruning
// argument for uniform-grid nearest-neighbour search, specialized to
// this project's AABB-overlap bin registration.
ClosestHit closestPointOnSoup(const TriangleAabbBins& bins, const Vec3& p);

// --- The smoothed signed distance field -------------------------------
//
// The layer march (`src/Layers.cpp`) steers off the EXACT distance-to-STL
// field, which has two properties that are today fought downstream, per
// point, by clamps and whole-stack condemnation: its gradient
// is discontinuous across the medial axis and every concave crease, and
// it retains every sub-cell facet/rib/sliver the grid cannot resolve.
// Smoothing moves that fight upstream: smooth the field once,
// at a radius tied to the local cell size, so the march sees a field
// smooth at the scale of the step it is taking.
//
// Definition: the smoothed
// field is the local average of the SIGNED distance over a ball of
// radius `r`, phi_r(p) = sum_k w_k * phi(p + r*o_k), sum_k w_k = 1, for a
// fixed rotation-fair quadrature {o_k, w_k} of the unit ball (see
// `sphere14Directions` below) evaluated on two radial shells (`r/2`, `r`) plus
// the centre point.
//
// Why THIS definition, and not a cheaper one:
// - NOT an exponential/power soft-min over the triangle list
//   (`-r log sum_i exp(-d_i/r)`): NOT tessellation-invariant -- a flat
//   wall split into N nearly-equidistant triangles yields a spurious
//   `-r log N` offset that grows with mesh density (fatal on real CAD).
//   Smoothing integrates over SPACE, not over the triangle list.
// - NOT morphological opening/closing via level-set offsetting: for an
//   exact SDF, `-r` then `+r` is the identity in convex regions; a
//   genuine closing needs the dilated set's own interior distance field,
//   more machinery than a convolution for the same effect.
// - The field convolved MUST be SIGNED. Convolving |phi| is wrong and
//   visibly destroys the geometry: |phi| is V-shaped across the surface,
//   so averaging lifts the wall's own minimum to ~r above zero and the
//   zero level set evaporates. (If an implementation appears to make the
//   surface vanish, this is the bug to suspect first.)
// - Rejected representations (voxel occupancy, grid-sampled narrow-band
//   SDF): both manufacture cell-scale staircase artefacts or need
//   pitch <= r/2 (~10x finer than the finest cell). This
//   query-time formulation has NO sampling pitch at all; `r` is chosen
//   on geometry alone.
//
// Sign source: the existing exact ray-parity test, `classifyVertices`
// (`src/CutData.hpp`), called ONCE per field evaluation over every
// stencil point of every query batched together -- never a second parity
// test (determinism and the grazing-ray tie-breaks stay in one place).
// vertexSolid[i] == true means INSIDE the solid (see classifyVertices'
// contract), so phi = (solid ? -1 : +1) * sqrt(closestPointOnSoup's
// distSq) at each stencil point.
struct SmoothedField {
    double radiusFrac = 1.0;  // r = radiusFrac * h_local (validated production default)
    bool singleShell = false; // cheaper one-shell quadrature; test-only, the mesher always uses two shells
};

// One direction of the degree-5 rotation-fair 14-point spherical rule: the
// cube's 6 face centres (weight 1/15) + 8 vertices (weight 3/40),
// verified to sum to 1 (6/15 + 8*3/40 = 0.4 + 0.6 = 1.0) by
// `test_smoothing`'s "quadrature sums to unity" check. Antipodally
// symmetric (every direction's negation is also in the set with the same
// weight) -- the property that makes a flat wall's smoothed value equal
// its exact value to round-off REGARDLESS of the chosen radial weights:
// for a locally linear field, opposite-direction pairs cancel the
// gradient term exactly, leaving only the field's own value at the
// centre. Do NOT replace with an axis-aligned-only 6-point stencil (face
// directions alone): that quadrature is anisotropic and makes the
// rounding direction-dependent.
struct SphereDir {
    Vec3 o;
    double w;
};
const std::vector<SphereDir>& sphere14Directions();

// Batched entry point (the primary API -- classification must be
// batched, see above): evaluates phi_r at every point in `p`, each with
// its own radius `r[i]` (a spatially varying radius, e.g. tied to a
// per-point local cell size, is supported, though it makes phi_r not a
// true distance function near refinement interfaces). `locationInMesh` is forwarded to `classifyVertices`
// unchanged.
// Optional `gradOut`: when non-null, filled with a gradient estimate of
// phi built from the SAME full-shell (`r`) stencil samples this call
// already evaluates -- no extra queries (the direction consumer
// needs `normalize(-grad phi_r)` without doubling the query cost).
//
// Derivation: for the 14-point rotation-fair rule above,
// `sum_k w_k * o_k (x) o_k == (1/3) I` exactly (a standard property of
// this degree-5 spherical design: it integrates x^2/y^2/z^2 to 1/3 and
// every cross term to 0, matching the unit sphere's own surface
// moments). So, to first order in `r`,
//   S := sum_k w_k * o_k * phi(p + r*o_k)
//      = phi(p) * sum_k w_k*o_k + r * sum_k w_k * o_k (o_k . grad phi)
//      = 0                      + r * (1/3) * grad phi
// (the first term vanishes because the direction set is antipodally
// symmetric with equal per-pair weight -- the same property that makes
// phi_r == phi on a flat wall regardless of radial weights). Hence
// `grad phi ~= 3 * S / r`, read directly off the full-shell samples
// already computed for phi_r itself.
std::vector<double> smoothedSignedDistanceBatch(const TriangleAabbBins& bins, const std::vector<Triangle>& tris,
                                                  const Vec3& locationInMesh, const std::vector<Vec3>& p,
                                                  const std::vector<double>& r, const SmoothedField& cfg,
                                                  std::vector<Vec3>* gradOut = nullptr);

// Scalar convenience form, expressed in terms of the batched API above
// (one code path) -- pays the batching overhead of a single-element
// vector; fine for tests and any non-hot-loop caller. The march
// must call the batched form directly over a whole front.
double smoothedSignedDistance(const TriangleAabbBins& bins, const std::vector<Triangle>& tris,
                                const Vec3& locationInMesh, const Vec3& p, double r, const SmoothedField& cfg);

} // namespace ninja
