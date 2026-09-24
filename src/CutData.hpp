// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "BaseMesh.hpp"
#include "Stl.hpp"

namespace ninja {

struct TriangleAabbBins; // Geometry.hpp

// --- small Vec3 arithmetic, shared by CutData and Cutter -----------------
inline Vec3 operator+(const Vec3& a, const Vec3& b) { return Vec3{a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return Vec3{a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(const Vec3& a, double s) { return Vec3{a.x * s, a.y * s, a.z * s}; }
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

// One STL entry of the `geometry` dict block (multiple allowed, one
// per named wall patch, in dict declaration order).
struct StlEntry {
    std::string stlPath; // resolved, relative to case dir
    std::string patchName;
    // Raw dict key (e.g. "box.stl"), before resolving to `stlPath` --
    // `type surface` refinement rules reference an STL by this raw
    // filename (the same key used as its `geometry` sub-dict name),
    // not by its output patch name.
    std::string rawFile;
};

// Dict-level geometry configuration: zero or more STL entries.
struct GeometryConfig {
    bool present = false;
    std::vector<StlEntry> stls; // declaration order == output patch order
    Vec3 locationInMesh;
    // `refinementGeometry`: STLs read ONLY by `type surface` refinement
    // rules -- never cut against, never a patch, never layered
    // (`patchName` unused).
    std::vector<StlEntry> sizing;
};

// Hash for a sorted pair of point indices (edge key).
struct EdgeKeyHash {
    std::size_t operator()(const std::pair<int, int>& k) const noexcept {
        return (static_cast<std::size_t>(static_cast<uint32_t>(k.first)) << 32) ^
               static_cast<std::size_t>(static_cast<uint32_t>(k.second));
    }
};

using EdgeKey = std::pair<int, int>;
inline EdgeKey makeEdgeKey(int a, int b) { return a < b ? EdgeKey{a, b} : EdgeKey{b, a}; }

// One shared per-edge cut intercept: its position plus the ordinal of
// the STL (`Triangle::solidId`) whose triangle produced the winning
// crossing. Classification runs on the combined triangle soup
//, but each intercept still remembers which solid
// it came from, for multi-STL patch assignment.
struct Intercept {
    Vec3 point;
    int solidId = 0;
};

// Shared, single-sourced cut geometry: one in/out status per grid vertex,
// one intercept per cut grid edge. Computed once; every cell consumes
// this data as a pure function.
struct CutData {
    std::vector<bool> vertexSolid; // indexed by base-mesh point index
    // Half-grid cut: the ON state, parallel to vertexSolid; empty in the
    // ordinary binary path. An ON vertex has vertexSolid == true, so
    // clipping and intercept placement treat it as "not a kept fluid corner"
    // (which is what stops a kept point and an intercept coinciding), but
    // CELL STATUS must NOT treat it as solid interior: doing so dilates the
    // removed region by the ON band and deletes whole cells at corners with
    // no cut face at all, gouging grid-aligned notches around convex
    // corners.
    std::vector<char> vertexOnSurface;
    std::unordered_map<EdgeKey, Intercept, EdgeKeyHash> edgeIntercept; // only cut edges
};

// Classifies every point in `points` in/out of the STL solid via
// ray-parity, referenced against `locationInMesh`'s parity so the
// caller doesn't need to know the STL's winding convention. `tris` may
// be the combined triangle soup of several STLs (union of solids,
// ) — classification does not need to know
// which solid each triangle belongs to. `prebuilt`, when given, must be
// buildTriangleAabbBins(tris) (default arguments) -- it saves rebuilding
// it on every call.
std::vector<bool> classifyVertices(const std::vector<Vec3>& points, const std::vector<Triangle>& tris,
                                    const Vec3& locationInMesh, const TriangleAabbBins* prebuilt = nullptr);

// Computes shared intercepts for every edge in `edges` (pairs of point
// indices into `points`) whose endpoint statuses differ.
std::unordered_map<EdgeKey, Intercept, EdgeKeyHash> computeEdgeIntercepts(
    const std::vector<Vec3>& points, const std::vector<Triangle>& tris,
    const std::vector<bool>& vertexSolid, const std::vector<EdgeKey>& edges);

// Collects the unique edges implied by a face->cells adjacency mesh's
// face point loops (consecutive pairs, deduplicated). General enough for
// any convex polyhedron (refined cells included).
std::vector<EdgeKey> collectEdges(const GeneratedMesh& mesh);

// --- The HALF-GRID offset cut (three-state, nearest snapping) -------------
//
// Builds CutData against the implicit offset solid
// { p : insideOriginal(p) OR (exists STL s: distToStl_s(p) <= stlThickness[s]) }
// instead of the plain STL solid. `tris` is the combined soup (used for
// parity classification and closestPointOnSoup's solidId attribution,
// exactly as the exact path uses it); `perStlTris`/`stlThickness` are
// indexed in parallel by STL declaration ordinal (== Triangle::solidId).
// A `stlThickness[s] == 0` STL contributes nothing beyond the plain parity
// term for that STL's own triangles, so an unlayered STL in a mixed case is
// cut at its true surface.
//
// Motivation: when t is a multiple
// of h the d = t level set lands exactly on grid planes, and `dist <= t` is
// then decided by floating-point noise -- e.g.
// |0.3-0.25| = 0.04999999999999999 (TRUE) vs |0.75-0.7| = 0.050000000000000044
// (FALSE), two ULPs apart in opposite directions, because 0.3 and 0.7 are not
// binary-exact. The offset region can then run one grid plane deeper on one
// side than the other, giving the two halves of a mirror-symmetric solid
// structurally different layer fronts.
//
// Fix: snap the surface to the nearest point of the HALF-GRID (grid vertices
// union edge midpoints) and classify each vertex into THREE states, with
// phi(v) = dist(v) - t and |grad phi| ~ 1:
//     |phi| < band  -> ON   (surface passes through this vertex)
//     phi   < 0     -> IN
//     phi   > 0     -> OUT
// `band` = h_local/4 is NOT an epsilon: it is the Voronoi radius of a
// half-grid node, the exact width at which "nearest" changes answer.
//
// The tie stops being load-bearing: rotcube's tied vertices sit at
// phi = -1.4e-17 and +4.4e-17, and BOTH land inside |phi| < 0.0125, so both
// are ON whichever side of zero the noise falls. Symmetry is preserved by
// construction rather than by a tie-break a representation change could flip.
// Quantization error also halves (reachable positions go from spacing h to
// h/2, worst-case displacement h/2 -> h/4).
//
// A+parity is untouched: the verdict is a pure function of position computed
// once per grid vertex and shared by every cell touching it, exactly as the
// binary parity is. Watertightness comes from the sharing, not the state count.
//
// Returns the binary "not kept" mask (ON and IN both count as solid, so an ON
// vertex is never a kept fluid corner -- which is what stops it being emitted
// as a kept point AND an intercept at the same position). `onOut` is sized to
// `points` and flags the ON state.
// Per-vertex thickness: where pointThickness[s] is non-empty, STL s
// offsets vertex i by pointThickness[s][i] (its local thickness field)
// instead of the constant stlThickness[s].
std::vector<bool> classifyVerticesOffsetHalfGrid(const std::vector<Vec3>& points,
                                                 const std::vector<Triangle>& tris,
                                                 const std::vector<std::vector<Triangle>>& perStlTris,
                                                 const Vec3& locationInMesh,
                                                 const std::vector<double>& stlThickness,
                                                 const std::vector<double>& bandPerPoint,
                                                 std::vector<char>& onOut,
                                                 const std::vector<std::vector<double>>& pointThickness = {});

// Half-grid cut statistics, reported by --cut-stats.
struct OffsetCutStats {
    int multiRootEdges = 0;
    // Edges cut at their midpoint (the nearest half-grid node).
    int quantizedIntercepts = 0;
    // Intercepts placed AT an ON grid vertex (vs at a midpoint).
    int onVertexIntercepts = 0;
};

std::unordered_map<EdgeKey, Intercept, EdgeKeyHash> computeEdgeInterceptsOffsetHalfGrid(
    const std::vector<Vec3>& points, const std::vector<Triangle>& tris, const std::vector<bool>& vertexSolid,
    const std::vector<char>& onSurface, const std::vector<EdgeKey>& edges, OffsetCutStats& stats);

} // namespace ninja
