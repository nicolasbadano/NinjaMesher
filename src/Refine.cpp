// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#include "Refine.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "Geometry.hpp"
#include "Par.hpp"

namespace ninja {

namespace {

// Global fine-grid index triple: a point is uniquely identified by its
// position in the finest grid implied by maxLevel, regardless of which
// leaf produced it (single-sourced, keyed by parent entity — here, the
// grid index itself). Base-grid vertices are the points with
// I%R==0 && J%R==0 && K%R==0.
struct Key {
    int i = 0;
    int j = 0;
    int k = 0;
    bool operator<(const Key& o) const {
        if (i != o.i) return i < o.i;
        if (j != o.j) return j < o.j;
        return k < o.k;
    }
    bool operator==(const Key& o) const { return i == o.i && j == o.j && k == o.k; }
};

int baseCellIndex(int i, int j, int k, int nx, int ny) { return i + j * nx + k * nx * ny; }

// The generated face loops
// (quads, later spliced into conforming polygons) are stored FLAT --
// one Key array + CSR offsets + parallel per-quad emittingCell/side
// arrays -- instead of a vector<RawQuad>{vector<Key>, int,
// std::string}: a measured 23.9 GB fishpassage peak sat inside
// assemble(serial), and the per-quad heap blocks and per-quad
// sideKey strings were a large slice of it. `side` is the domain-side
// ordinal 0..5 (index into sideKeys6) or -1 for a non-boundary quad --
// the former sideKey string was always one of those six literals.
struct QuadStore {
    std::vector<Key> keys;      // flat loops
    std::vector<int> offsets{0};
    std::vector<int> cell;      // emitting cell per quad
    std::vector<signed char> side; // 0..5 domain side, -1 internal
    int size() const { return static_cast<int>(cell.size()); }
    int count(int q) const {
        return offsets[static_cast<std::size_t>(q) + 1] - offsets[static_cast<std::size_t>(q)];
    }
    const Key* ptsOf(int q) const { return keys.data() + offsets[static_cast<std::size_t>(q)]; }
    void append(const Key* pts, int n, int emittingCell, int sideOrd) {
        keys.insert(keys.end(), pts, pts + n);
        offsets.push_back(checkedInt(keys.size()));
        cell.push_back(emittingCell);
        side.push_back(static_cast<signed char>(sideOrd));
    }
};

int trailingZeros(int v) {
    int c = 0;
    while ((v & 1) == 0) {
        v >>= 1;
        ++c;
    }
    return c;
}

// One octree leaf: `i0,j0,k0` is its minimum corner in GLOBAL fine-grid
// units (the same convention as Key); its extent is `size = R >> level`
// fine cells per axis; `level` is the leaf's depth (0 == the whole base
// cell, unsplit).
struct Leaf {
    int i0 = 0, j0 = 0, k0 = 0;
    int level = 0;
};

// The 4 corner keys of sub-quad (p,q) of a leaf face -- the EXACT
// corner geometry the serial assembly's emitDir uses,
// shared by assembleRefined's quad
// emission and refineLocal's block-local corner/edge enumeration so
// the two can never drift.
std::array<Key, 4> subQuadCorners(int dir, int I0, int J0, int K0, int step, int F, int p, int q) {
    const int sub = step / F;
    std::array<Key, 4> corners{};
    switch (dir) {
        case 0: { // -x, fixed I=I0, (p=y,q=z)
            const int j0 = J0 + p * sub, j1 = j0 + sub, k0 = K0 + q * sub, k1 = k0 + sub;
            corners = {Key{I0, j0, k0}, Key{I0, j0, k1}, Key{I0, j1, k1}, Key{I0, j1, k0}};
            break;
        }
        case 1: { // +x, fixed I=I0+step
            const int I1 = I0 + step;
            const int j0 = J0 + p * sub, j1 = j0 + sub, k0 = K0 + q * sub, k1 = k0 + sub;
            corners = {Key{I1, j0, k0}, Key{I1, j1, k0}, Key{I1, j1, k1}, Key{I1, j0, k1}};
            break;
        }
        case 2: { // -y, fixed J=J0, (p=x,q=z)
            const int i0 = I0 + p * sub, i1 = i0 + sub, k0 = K0 + q * sub, k1 = k0 + sub;
            corners = {Key{i0, J0, k0}, Key{i1, J0, k0}, Key{i1, J0, k1}, Key{i0, J0, k1}};
            break;
        }
        case 3: { // +y, fixed J=J0+step
            const int J1 = J0 + step;
            const int i0 = I0 + p * sub, i1 = i0 + sub, k0 = K0 + q * sub, k1 = k0 + sub;
            corners = {Key{i0, J1, k0}, Key{i0, J1, k1}, Key{i1, J1, k1}, Key{i1, J1, k0}};
            break;
        }
        case 4: { // -z, fixed K=K0, (p=x,q=y)
            const int i0 = I0 + p * sub, i1 = i0 + sub, j0 = J0 + q * sub, j1 = j0 + sub;
            corners = {Key{i0, j0, K0}, Key{i0, j1, K0}, Key{i1, j1, K0}, Key{i1, j0, K0}};
            break;
        }
        case 5: { // +z, fixed K=K0+step
            const int K1 = K0 + step;
            const int i0 = I0 + p * sub, i1 = i0 + sub, j0 = J0 + q * sub, j1 = j0 + sub;
            corners = {Key{i0, j0, K1}, Key{i1, j0, K1}, Key{i1, j1, K1}, Key{i0, j1, K1}};
            break;
        }
        default:
            break;
    }
    return corners;
}

// The corner keys of a quad set, one sorted copy per axis ordered so
// the keys on any grid line along that axis are contiguous and
// ascending along it: the conformance splice then reads the corners
// strictly inside a quad edge as one range, instead of probing every
// fine lattice point of the edge (R of them, 2048 at maxLevel 11).
class CornerLines {
public:
    explicit CornerLines(std::vector<Key> keys) {
        std::sort(keys.begin(), keys.end());
        keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
        for (int a = 0; a < 3; ++a) {
            byAxis_[a] = keys;
            std::sort(byAxis_[a].begin(), byAxis_[a].end(),
                      [a](const Key& u, const Key& v) { return lineOrder(a, u) < lineOrder(a, v); });
        }
    }

    // Appends the corners strictly between the ends of the axis-aligned
    // segment a->b, in order from a to b.
    void appendInterior(const Key& a, const Key& b, std::vector<Key>& out) const {
        const int ax = a.i != b.i ? 0 : (a.j != b.j ? 1 : (a.k != b.k ? 2 : -1));
        if (ax < 0) return;
        const std::array<int, 3> la = lineOrder(ax, a), lb = lineOrder(ax, b);
        std::array<int, 3> lo = la, hi = la;
        lo[2] = std::min(la[2], lb[2]) + 1;
        hi[2] = std::max(la[2], lb[2]) - 1;
        if (lo[2] > hi[2]) return;
        const std::vector<Key>& v = byAxis_[ax];
        auto cmpLo = [ax](const Key& k, const std::array<int, 3>& t) { return lineOrder(ax, k) < t; };
        auto cmpHi = [ax](const std::array<int, 3>& t, const Key& k) { return t < lineOrder(ax, k); };
        const auto first = std::lower_bound(v.begin(), v.end(), lo, cmpLo);
        const auto last = std::upper_bound(first, v.end(), hi, cmpHi);
        if (la[2] < lb[2]) {
            out.insert(out.end(), first, last);
        } else {
            out.insert(out.end(), std::make_reverse_iterator(last), std::make_reverse_iterator(first));
        }
    }

private:
    // (the two fixed coordinates, then the coordinate along `axis`)
    static std::array<int, 3> lineOrder(int axis, const Key& k) {
        if (axis == 0) return {k.j, k.k, k.i};
        if (axis == 1) return {k.i, k.k, k.j};
        return {k.i, k.j, k.k};
    }
    std::array<std::vector<Key>, 3> byAxis_;
};

// --- Static 3D block decomposition of the base grid ---------------------
// (parDims3 = MPI_Dims_create.) Rank layout
// convention: rank = cx + dims[0]*(cy + dims[1]*cz). Empty blocks (np
// larger than an axis' cell count) sit at the HIGH end of each axis by
// parBlockRange's construction, so an empty neighbour block implies the
// domain effectively ends at this rank's boundary.
struct Decomp {
    int n[3] = {0, 0, 0};     // base cells per axis
    int dims[3] = {1, 1, 1};  // process grid
    int coord[3] = {0, 0, 0}; // my coords in it
    int ownLo[3] = {0, 0, 0}; // my base-cell block [lo, hi)
    int ownHi[3] = {0, 0, 0};
    int extLo[3] = {0, 0, 0}; // block + one-base-cell halo, clamped
    int extHi[3] = {0, 0, 0};

    bool ownEmpty() const {
        return ownHi[0] <= ownLo[0] || ownHi[1] <= ownLo[1] || ownHi[2] <= ownLo[2];
    }
    int rankAt(int cx, int cy, int cz) const { return cx + dims[0] * (cy + dims[1] * cz); }
};

Decomp makeDecomp(int nx, int ny, int nz) {
    Decomp d;
    d.n[0] = nx;
    d.n[1] = ny;
    d.n[2] = nz;
    const std::array<int, 3> dims = parDims3();
    d.dims[0] = dims[0];
    d.dims[1] = dims[1];
    d.dims[2] = dims[2];
    const int rank = parRank();
    d.coord[0] = rank % d.dims[0];
    d.coord[1] = (rank / d.dims[0]) % d.dims[1];
    d.coord[2] = rank / (d.dims[0] * d.dims[1]);
    for (int a = 0; a < 3; ++a) {
        std::size_t lo = 0, hi = 0;
        parBlockRange(static_cast<std::size_t>(d.n[a]), d.coord[a], d.dims[a], lo, hi);
        d.ownLo[a] = static_cast<int>(lo);
        d.ownHi[a] = static_cast<int>(hi);
        d.extLo[a] = std::max(0, d.ownLo[a] - 1);
        d.extHi[a] = std::min(d.n[a], d.ownHi[a] + 1);
    }
    if (d.ownEmpty()) {
        for (int a = 0; a < 3; ++a) {
            d.extLo[a] = d.ownLo[a];
            d.extHi[a] = d.ownLo[a];
        }
    }
    return d;
}

// Block (process-grid) index owning base cell `c` along an axis of `n`
// cells split over `P` blocks -- the exact inverse of parBlockRange.
int blockOfCell(int c, int n, int P) {
    const int base = n / P;
    const int rem = n % P;
    if (base == 0) return c; // c < rem by construction
    const int cut = (base + 1) * rem;
    return c < cut ? c / (base + 1) : rem + (c - cut) / base;
}

// Per-octant octree over a box of base cells (a rank's extended =
// owned + halo box, or the whole domain for the rank-0 assembly).
// Leaf coordinates are GLOBAL fine-grid units throughout. The fine-cell
// lookup walks down from the base cell's root (split children are
// created consecutively, octant = dx + 2*dy + 4*dz), so memory scales
// with the number of leaves, not with the fine-grid volume (R^3 per
// base cell, which at maxLevel 11 is 8.6e9).
class LocalOctree {
public:
    LocalOctree(const int baseLo[3], const int baseHi[3], int R) : R_(R) {
        for (int a = 0; a < 3; ++a) {
            bLo_[a] = baseLo[a];
            bDim_[a] = std::max(0, baseHi[a] - baseLo[a]);
            loF_[a] = baseLo[a] * R;
            hiF_[a] = baseHi[a] * R;
        }
        roots_.assign(static_cast<std::size_t>(bDim_[0]) * static_cast<std::size_t>(bDim_[1]) *
                          static_cast<std::size_t>(bDim_[2]),
                      -1);
        // One initial leaf per base cell, level 0.
        for (int bk = 0; bk < bDim_[2]; ++bk) {
            for (int bj = 0; bj < bDim_[1]; ++bj) {
                for (int bi = 0; bi < bDim_[0]; ++bi) {
                    roots_[rootFlat(bi, bj, bk)] =
                        newLeaf((bLo_[0] + bi) * R, (bLo_[1] + bj) * R, (bLo_[2] + bk) * R, 0);
                }
            }
        }
    }

    int R() const { return R_; }
    int leafCount() const { return static_cast<int>(leaves_.size()); }

    std::vector<int> liveIds() const {
        std::vector<int> out;
        out.reserve(leaves_.size());
        for (std::size_t i = 0; i < leaves_.size(); ++i) {
            if (live_[i]) out.push_back(static_cast<int>(i));
        }
        return out;
    }

    const Leaf& leaf(int id) const { return leaves_[static_cast<std::size_t>(id)]; }
    int leafSize(int id) const { return R_ >> leaves_[static_cast<std::size_t>(id)].level; }

    // Leaf id covering GLOBAL fine cell (i,j,k); -1 outside the box
    // (which subsumes "outside the domain" -- both mean "no neighbour
    // visible", exactly the serial semantics at the domain boundary;
    // owned leaves only ever probe one fine cell out, which is inside
    // the extended box unless it is outside the domain).
    int at(int i, int j, int k) const {
        if (i < loF_[0] || i >= hiF_[0] || j < loF_[1] || j >= hiF_[1] || k < loF_[2] || k >= hiF_[2]) {
            return -1;
        }
        int id = roots_[rootFlat(i / R_ - bLo_[0], j / R_ - bLo_[1], k / R_ - bLo_[2])];
        while (firstChild_[static_cast<std::size_t>(id)] >= 0) {
            const Leaf& L = leaves_[static_cast<std::size_t>(id)];
            const int half = (R_ >> L.level) >> 1;
            id = firstChild_[static_cast<std::size_t>(id)] + (i - L.i0 >= half ? 1 : 0) +
                 (j - L.j0 >= half ? 2 : 0) + (k - L.k0 >= half ? 4 : 0);
        }
        return id;
    }

    // Max level of the leaves covering any fine cell of the half-open
    // GLOBAL fine box [lo, hi) clipped to this tree's box; -1 when the
    // clipped box is empty. Equal to the max of at() over every fine
    // cell of the box, without visiting them one by one.
    int maxLevelIn(const int lo[3], const int hi[3]) const {
        int c[3], d[3];
        for (int a = 0; a < 3; ++a) {
            c[a] = std::max(lo[a], loF_[a]);
            d[a] = std::min(hi[a], hiF_[a]);
            if (c[a] >= d[a]) return -1;
        }
        int best = -1;
        std::vector<int> stack;
        for (int bk = c[2] / R_; bk <= (d[2] - 1) / R_; ++bk) {
            for (int bj = c[1] / R_; bj <= (d[1] - 1) / R_; ++bj) {
                for (int bi = c[0] / R_; bi <= (d[0] - 1) / R_; ++bi) {
                    stack.push_back(roots_[rootFlat(bi - bLo_[0], bj - bLo_[1], bk - bLo_[2])]);
                }
            }
        }
        while (!stack.empty()) {
            const int id = stack.back();
            stack.pop_back();
            const int fc = firstChild_[static_cast<std::size_t>(id)];
            if (fc < 0) {
                best = std::max(best, leaves_[static_cast<std::size_t>(id)].level);
                continue;
            }
            for (int o = 0; o < 8; ++o) {
                const Leaf& C = leaves_[static_cast<std::size_t>(fc + o)];
                const int s = R_ >> C.level;
                if (C.i0 < d[0] && C.i0 + s > c[0] && C.j0 < d[1] && C.j0 + s > c[1] && C.k0 < d[2] &&
                    C.k0 + s > c[2]) {
                    stack.push_back(fc + o);
                }
            }
        }
        return best;
    }

    // Live leaves inside the base-cell-aligned GLOBAL fine box [lo, hi),
    // in depth-first order from the base cells taken k-major.
    void leavesIn(const int lo[3], const int hi[3], std::vector<LeafRec>& out) const {
        std::vector<int> stack;
        for (int bk = lo[2] / R_; bk < hi[2] / R_; ++bk) {
            for (int bj = lo[1] / R_; bj < hi[1] / R_; ++bj) {
                for (int bi = lo[0] / R_; bi < hi[0] / R_; ++bi) {
                    stack.push_back(roots_[rootFlat(bi - bLo_[0], bj - bLo_[1], bk - bLo_[2])]);
                    while (!stack.empty()) {
                        const int id = stack.back();
                        stack.pop_back();
                        const int fc = firstChild_[static_cast<std::size_t>(id)];
                        if (fc < 0) {
                            const Leaf& L = leaves_[static_cast<std::size_t>(id)];
                            out.push_back(LeafRec{L.i0, L.j0, L.k0, L.level});
                            continue;
                        }
                        for (int o = 7; o >= 0; --o) stack.push_back(fc + o);
                    }
                }
            }
        }
    }

    // Splits the leaf containing fine cell (i,j,k) until it is at least
    // `level` deep.
    void refineTo(int i, int j, int k, int level) {
        for (int id = at(i, j, k); leaves_[static_cast<std::size_t>(id)].level < level; id = at(i, j, k)) {
            split(id);
        }
    }

    // Splits leaf `id` into 8 children (level+1).
    void split(int id) {
        const Leaf parent = leaves_[static_cast<std::size_t>(id)];
        live_[static_cast<std::size_t>(id)] = false;
        const int childSize = R_ >> (parent.level + 1);
        int first = -1;
        for (int dz = 0; dz <= 1; ++dz) {
            for (int dy = 0; dy <= 1; ++dy) {
                for (int dx = 0; dx <= 1; ++dx) {
                    const int cid = newLeaf(parent.i0 + dx * childSize, parent.j0 + dy * childSize,
                                            parent.k0 + dz * childSize, parent.level + 1);
                    if (first < 0) first = cid;
                }
            }
        }
        firstChild_[static_cast<std::size_t>(id)] = first;
    }

private:
    std::size_t rootFlat(int bi, int bj, int bk) const {
        return static_cast<std::size_t>(bi) + static_cast<std::size_t>(bj) * static_cast<std::size_t>(bDim_[0]) +
               static_cast<std::size_t>(bk) * static_cast<std::size_t>(bDim_[0]) *
                   static_cast<std::size_t>(bDim_[1]);
    }

    int newLeaf(int i0, int j0, int k0, int level) {
        leaves_.push_back(Leaf{i0, j0, k0, level});
        live_.push_back(true);
        firstChild_.push_back(-1);
        return static_cast<int>(leaves_.size()) - 1;
    }

    int R_;
    int bLo_[3] = {0, 0, 0}, bDim_[3] = {0, 0, 0};
    int loF_[3] = {0, 0, 0}, hiF_[3] = {0, 0, 0};
    std::vector<Leaf> leaves_;
    std::vector<char> live_;
    std::vector<int> firstChild_; // -1 for a live leaf
    std::vector<int> roots_;      // per base cell of the box
};

// --- The 2:1-grading halo exchange.
// Exchanges the LEAVES of the one-base-cell boundary layer with
// face-neighbour ranks, one axis phase at a time (x, then y, then z),
// each phase's slab spanning the EXTENDED range of the already-exchanged
// axes -- the standard trick that delivers edge/corner-diagonal halo
// cells through face-neighbour messages only. Receivers IMPOSE the
// levels by splitting halo leaves until each received leaf's corner is
// at least that deep; that covers the whole received leaf, since a
// shallower local leaf overlapping it would contain it, corner
// included (levels only ever grow, so imposition is monotone and
// exact). The conformance-splice corner keys need no separate message
// type: they are a pure function of these same exchanged halo levels
// (Refine.hpp).
void exchangeHalos(LocalOctree& tree, const Decomp& dc, int R) {
    if (dc.ownEmpty()) {
        return; // no cells, no partners (see Decomp: empty blocks stack at the high end)
    }
    std::vector<LeafRec> sendBuf, recvBuf;
    for (int axis = 0; axis < 3; ++axis) {
        // Tangential fine ranges for this phase: already-exchanged axes
        // (a < axis) use the EXTENDED range, later axes the OWN range.
        int tLo[3], tHi[3];
        for (int a = 0; a < 3; ++a) {
            if (a < axis) {
                tLo[a] = dc.extLo[a] * R;
                tHi[a] = dc.extHi[a] * R;
            } else {
                tLo[a] = dc.ownLo[a] * R;
                tHi[a] = dc.ownHi[a] * R;
            }
        }

        // Partner existence: neighbour coord in range AND its block
        // non-empty (empty blocks stack at the axis' high end, so an
        // empty +side neighbour implies ownHi == n: a domain boundary).
        auto partnerRank = [&](int delta) {
            const int c = dc.coord[axis] + delta;
            if (c < 0 || c >= dc.dims[axis]) return -1;
            std::size_t lo = 0, hi = 0;
            parBlockRange(static_cast<std::size_t>(dc.n[axis]), c, dc.dims[axis], lo, hi);
            if (hi <= lo) return -1;
            int coords[3] = {dc.coord[0], dc.coord[1], dc.coord[2]};
            coords[axis] = c;
            return dc.rankAt(coords[0], coords[1], coords[2]);
        };
        const int prev = partnerRank(-1);
        const int next = partnerRank(+1);

        auto box = [&](int fineLo, int* lo, int* hi) {
            for (int a = 0; a < 3; ++a) {
                lo[a] = tLo[a];
                hi[a] = tHi[a];
            }
            lo[axis] = fineLo;
            hi[axis] = fineLo + R;
        };
        // Symmetric Sendrecv of a leaf list: counts first, then records
        // (deadlock-free by construction).
        auto shift = [&](int dest, int src, int tag) {
            std::uint64_t nSend = sendBuf.size(), nRecv = 0;
            parSendRecvBytes(dest, &nSend, dest >= 0 ? sizeof(nSend) : 0, src, &nRecv,
                             src >= 0 ? sizeof(nRecv) : 0, tag);
            recvBuf.assign(static_cast<std::size_t>(nRecv), LeafRec{});
            parSendRecvBytes(dest, sendBuf.data(), dest >= 0 ? sendBuf.size() * sizeof(LeafRec) : 0, src,
                             recvBuf.data(), recvBuf.size() * sizeof(LeafRec), tag + 1000);
            for (const LeafRec& L : recvBuf) tree.refineTo(L.i0, L.j0, L.k0, L.level);
        };
        int sLo[3], sHi[3];

        // Shift in +axis direction: everyone sends its TOP owned layer
        // to `next` and receives `prev`'s top layer into its -side halo.
        sendBuf.clear();
        box(dc.ownHi[axis] * R - R, sLo, sHi);
        if (next >= 0) tree.leavesIn(sLo, sHi, sendBuf);
        shift(next, prev, 100 + axis);

        // Shift in -axis direction: send my BOTTOM owned layer to
        // `prev`, receive `next`'s bottom layer into my +side halo.
        sendBuf.clear();
        box(dc.ownLo[axis] * R, sLo, sHi);
        if (prev >= 0) tree.leavesIn(sLo, sHi, sendBuf);
        shift(prev, next, 200 + axis);
    }
}

} // namespace

LocalRefine refineLocal(const MeshConfig& cfg, const std::vector<RefineRegion>& regions,
                        const std::unordered_map<std::string, std::vector<Triangle>>& stlTriangles,
                        bool needCutEntities) {
    const int nx = cfg.nx;
    const int ny = cfg.ny;
    const int nz = cfg.nz;
    if (nx <= 0 || ny <= 0 || nz <= 0) {
        throw std::runtime_error("Refine error: n must be positive in all directions");
    }

    // maxLevel is a pure function of the requested region levels.
    int maxLevel = 0;
    for (const RefineRegion& r : regions) {
        maxLevel = std::max(maxLevel, r.level);
    }
    const int R = 1 << maxLevel;

    // Distance-mode regions: bin each named STL's triangles (replicated
    // STLs, the decided model -- every rank holds the full soup).
    std::vector<TriangleAabbBins> surfaceBins(regions.size());
    for (std::size_t r = 0; r < regions.size(); ++r) {
        if (regions[r].type != "surface") continue;
        auto it = stlTriangles.find(regions[r].stlKey);
        if (it == stlTriangles.end()) {
            throw std::runtime_error("Refine error: surface region references unknown STL key '" +
                                      regions[r].stlKey + "' (not declared in 'geometry')");
        }
        surfaceBins[r] = buildTriangleAabbBins(it->second);
    }

    const Decomp dc = makeDecomp(nx, ny, nz);
    LocalOctree tree(dc.extLo, dc.extHi, R);

    auto centroidOf = [&](int i0, int j0, int k0, int size) {
        return Vec3{cfg.min.x + (cfg.max.x - cfg.min.x) * (i0 + size * 0.5) / (nx * R),
                    cfg.min.y + (cfg.max.y - cfg.min.y) * (j0 + size * 0.5) / (ny * R),
                    cfg.min.z + (cfg.max.z - cfg.min.z) * (k0 + size * 0.5) / (nz * R)};
    };

    for (const RefineRegion& reg : regions) {
        if (reg.type != "box" && reg.type != "sphere" && reg.type != "surface") {
            throw std::runtime_error("Refine error: unknown region type '" + reg.type + "'");
        }
    }
    // True when some region deeper than `level` contains the leaf's
    // centroid, i.e. when the leaf must split. Cheap region types are
    // tested first; the first hit decides.
    auto wantsSplit = [&](int i0, int j0, int k0, int size, int level) {
        const Vec3 c = centroidOf(i0, j0, k0, size);
        for (int pass = 0; pass < 2; ++pass) {
            for (std::size_t r = 0; r < regions.size(); ++r) {
                const RefineRegion& reg = regions[r];
                if (reg.level <= level || (reg.type == "surface") != (pass == 1)) continue;
                bool inside = false;
                if (reg.type == "box") {
                    inside = c.x >= reg.min.x && c.x <= reg.max.x && c.y >= reg.min.y && c.y <= reg.max.y &&
                             c.z >= reg.min.z && c.z <= reg.max.z;
                } else if (reg.type == "sphere") {
                    const Vec3 d = c - reg.centre;
                    inside = dot(d, d) <= reg.radius * reg.radius;
                } else {
                    inside = anyTriangleWithin(surfaceBins[r], c, reg.distance);
                }
                if (inside) return true;
            }
        }
        return false;
    };

    auto ownsLeaf = [&](const Leaf& L) {
        const int bi = L.i0 / R;
        const int bj = L.j0 / R;
        const int bk = L.k0 / R;
        return bi >= dc.ownLo[0] && bi < dc.ownHi[0] && bj >= dc.ownLo[1] && bj < dc.ownHi[1] &&
               bk >= dc.ownLo[2] && bk < dc.ownHi[2];
    };

    // --- Iterative mark -> split -> grade -> repeat, to GLOBAL
    // stability: the serial loop restricted to owned leaves,
    // with the halo exchange at the top of every grading iteration and
    // Allreduce fixpoint termination. At
    // np=1 every collective is a no-op and this is exactly the serial
    // loop. Multi-rank: both marking (once per owned leaf per outer
    // pass) and grading (a monotone split-only closure) are pure, so
    // per-rank chaotic iteration converges to the same unique least
    // fixpoint as the serial sweep -- the final leaf set is
    // partition-invariant by construction (and gated byte-identically).
    std::vector<char> settled; // per leaf id: Pass A found no split needed
    // NINJA_STAGE_TIMES: one line per marking/grading sweep.
    static const bool showProgress = std::getenv("NINJA_STAGE_TIMES") != nullptr;
    const auto t0 = std::chrono::steady_clock::now();
    auto progress = [&](const char* what, std::size_t live, std::size_t nSplit) {
        if (!showProgress || parRank() != 0) return;
        const double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::cout << "refine " << what << " rank 0: " << live << " live leaves, " << nSplit << " split, " << t
                  << " s\n" << std::flush;
    };
    bool outerChanged = true;
    while (outerChanged) {
        int changed = 0;

        // Pass A: region/distance marking (owned leaves only -- each
        // leaf has exactly one owner, so the -- possibly expensive --
        // distance predicate is evaluated exactly once globally; a leaf
        // found not to need a split keeps that answer, the predicate
        // being a pure function of the leaf).
        {
            const std::vector<int> ids = tree.liveIds();
            std::vector<char> split(ids.size(), 0);
            settled.resize(std::max(settled.size(), static_cast<std::size_t>(tree.leafCount())), 0);
#pragma omp parallel for schedule(dynamic, 256)
            for (std::size_t n = 0; n < ids.size(); ++n) {
                const int id = ids[n];
                if (settled[static_cast<std::size_t>(id)] != 0) continue;
                const Leaf& L = tree.leaf(id);
                if (!ownsLeaf(L)) continue;
                split[n] = wantsSplit(L.i0, L.j0, L.k0, tree.leafSize(id), L.level) ? 1 : 0;
                if (split[n] == 0) settled[static_cast<std::size_t>(id)] = 1;
            }
            std::size_t nSplit = 0;
            for (std::size_t n = 0; n < ids.size(); ++n) {
                if (split[n] == 0) continue;
                tree.split(ids[n]);
                changed = 1;
                ++nSplit;
            }
            progress("mark", ids.size(), nSplit);
        }

        // Pass B: 2:1 grading on owned leaves, halo-exchanged and
        // iterated to a GLOBAL fixpoint.
        for (;;) {
            exchangeHalos(tree, dc, R);
            const std::vector<int> ids = tree.liveIds();
            std::vector<char> split(ids.size(), 0);
#pragma omp parallel for schedule(dynamic, 512)
            for (std::size_t n = 0; n < ids.size(); ++n) {
                const int id = ids[n];
                const Leaf& L = tree.leaf(id);
                if (!ownsLeaf(L)) continue;
                const int size = tree.leafSize(id);
                // Max level over the one-fine-cell slab beyond each face
                // (grading only needs the MAX level touching the face).
                int maxNeighLevel = L.level;
                for (int dir = 0; dir < 6; ++dir) {
                    int lo[3] = {L.i0, L.j0, L.k0};
                    int hi[3] = {L.i0 + size, L.j0 + size, L.k0 + size};
                    const int a = dir / 2;
                    if (dir % 2 == 0) {
                        hi[a] = lo[a];
                        lo[a] -= 1;
                    } else {
                        lo[a] = hi[a];
                        hi[a] += 1;
                    }
                    maxNeighLevel = std::max(maxNeighLevel, tree.maxLevelIn(lo, hi));
                }
                if (maxNeighLevel - 1 > L.level) split[n] = 1;
            }
            std::vector<int> toSplit;
            for (std::size_t n = 0; n < ids.size(); ++n) {
                if (split[n] != 0) toSplit.push_back(ids[n]);
            }
            for (int id : toSplit) {
                tree.split(id);
            }
            progress("grade", ids.size(), toSplit.size());
            int g = toSplit.empty() ? 0 : 1;
            changed |= g;
            parAllreduceSumInts(&g, 1);
            if (g == 0) break;
        }

        int oc = changed;
        parAllreduceSumInts(&oc, 1);
        outerChanged = oc != 0;
    }

    LocalRefine out;
    out.R = R;
    for (int id : tree.liveIds()) {
        const Leaf& L = tree.leaf(id);
        if (ownsLeaf(L)) {
            out.ownedLeaves.push_back(LeafRec{L.i0, L.j0, L.k0, L.level});
        }
    }
    if (!needCutEntities || dc.ownEmpty()) {
        return out;
    }

    // --- Block-local cut-entity enumeration ------------------------------
    // Reproduces, for this rank's owned faces, exactly the spliced
    // point loops the serial assembly emits -- so the union over ranks
    // of {points} covers every serial mesh point, and the disjoint
    // union of {owned edges} IS the serial `collectEdges` edge set
    // (each face is emitted by its leaves' owner ranks; both sides of a
    // shared face carry the same spliced point set, hence the same
    // consecutive sub-segments).
    const int nxF = nx * R;
    const int nyF = ny * R;
    // Global lattice flat of a point key: ascending flat == ascending
    // (k, j, i), the point order.
    auto flatOf = [&](const Key& kk) { return fineLatticeFlat(kk.i, kk.j, kk.k, nxF, nyF); };

    // Face subdivision factor F for a leaf face, EXACTLY the serial
    // assembly's quadrant-probe logic (a probe outside the extended box
    // reads -1 == "no neighbour", which for OWNED leaves only happens
    // outside the DOMAIN -- serial's own convention; for halo leaves'
    // outward faces it may under-subdivide, harmless: those faces'
    // corners cannot lie on any owned face's edges, being two planes
    // away from the owned box).
    auto faceInfo = [&](const Leaf& L, int dir, bool& isBoundary) {
        const int size = R >> L.level;
        int probeI = L.i0, probeJ = L.j0, probeK = L.k0;
        isBoundary = false;
        switch (dir) {
            case 0: isBoundary = (L.i0 == 0); probeI = L.i0 - 1; break;
            case 1: isBoundary = (L.i0 + size == nx * R); probeI = L.i0 + size; break;
            case 2: isBoundary = (L.j0 == 0); probeJ = L.j0 - 1; break;
            case 3: isBoundary = (L.j0 + size == ny * R); probeJ = L.j0 + size; break;
            case 4: isBoundary = (L.k0 == 0); probeK = L.k0 - 1; break;
            case 5: isBoundary = (L.k0 + size == nz * R); probeK = L.k0 + size; break;
            default: break;
        }
        if (isBoundary) return 1;
        int F = 1;
        if (size >= 2) {
            const int half = size / 2;
            int ids[4];
            int n = 0;
            for (int a = 0; a <= half; a += half) {
                for (int b = 0; b <= half; b += half) {
                    int qi = probeI, qj = probeJ, qk = probeK;
                    switch (dir) {
                        case 0:
                        case 1:
                            qj = L.j0 + a;
                            qk = L.k0 + b;
                            break;
                        case 2:
                        case 3:
                            qi = L.i0 + a;
                            qk = L.k0 + b;
                            break;
                        default:
                            qi = L.i0 + a;
                            qj = L.j0 + b;
                            break;
                    }
                    ids[n++] = tree.at(qi, qj, qk);
                }
            }
            const bool allSame = ids[0] == ids[1] && ids[0] == ids[2] && ids[0] == ids[3];
            if (!allSame) F = 2;
        }
        return F;
    };

    // Pass 1: corner bitmap over ALL extended-box leaves (owned + halo)
    // -- this realizes the corner-key exchange: neighbour
    // ranks' quad corners are REGENERATED here from the exchanged halo
    // levels (a pure function of leaf geometry) instead of being
    // messaged as explicit keys.
    std::vector<Key> cornerKeys;
    const std::vector<int> liveAll = tree.liveIds();
    for (int id : liveAll) {
        const Leaf& L = tree.leaf(id);
        const int size = R >> L.level;
        for (int dir = 0; dir < 6; ++dir) {
            bool isBoundary = false;
            const int F = faceInfo(L, dir, isBoundary);
            for (int p = 0; p < F; ++p) {
                for (int q = 0; q < F; ++q) {
                    const std::array<Key, 4> c = subQuadCorners(dir, L.i0, L.j0, L.k0, size, F, p, q);
                    cornerKeys.insert(cornerKeys.end(), c.begin(), c.end());
                }
            }
        }
    }
    const CornerLines corners(std::move(cornerKeys));

    // Pass 2: owned leaves' quads, spliced against the corner set
    // (the serial conformance rule verbatim: insert every fine-grid
    // point strictly inside the edge that is used as a corner anywhere,
    // sorted along the edge); consecutive loop pairs are the serial
    // edge sub-segments. Edge OWNERSHIP (each global edge computed by
    // exactly one rank, so per-edge work/stats partition exactly): the
    // rank owning the base cell at the edge's min corner, +side on
    // shared planes, clamped at the domain max plane.
    auto edgeOwnedByMe = [&](const Key& u, const Key& v) {
        for (int a = 0; a < 3; ++a) {
            const int ua = a == 0 ? u.i : (a == 1 ? u.j : u.k);
            const int va = a == 0 ? v.i : (a == 1 ? v.j : v.k);
            int cellA;
            if (ua == va) {
                cellA = std::min(dc.n[a] - 1, ua / R);
            } else {
                cellA = std::min(ua, va) / R;
            }
            if (blockOfCell(cellA, dc.n[a], dc.dims[a]) != dc.coord[a]) {
                return false;
            }
        }
        return true;
    };

    std::vector<Key> pointKeys;
    std::vector<std::pair<std::int64_t, std::int64_t>> edgePairs;
    std::vector<Key> loop;
    for (int id : liveAll) {
        const Leaf& L = tree.leaf(id);
        if (!ownsLeaf(L)) continue;
        const int size = R >> L.level;
        for (int dir = 0; dir < 6; ++dir) {
            bool isBoundary = false;
            const int F = faceInfo(L, dir, isBoundary);
            for (int p = 0; p < F; ++p) {
                for (int q = 0; q < F; ++q) {
                    const std::array<Key, 4> c = subQuadCorners(dir, L.i0, L.j0, L.k0, size, F, p, q);
                    loop.clear();
                    for (std::size_t e = 0; e < 4; ++e) {
                        loop.push_back(c[e]);
                        corners.appendInterior(c[e], c[(e + 1) % 4], loop);
                    }
                    const std::size_t n = loop.size();
                    for (std::size_t e = 0; e < n; ++e) {
                        const Key& u = loop[e];
                        const Key& v = loop[(e + 1) % n];
                        pointKeys.push_back(u);
                        if (edgeOwnedByMe(u, v)) {
                            std::int64_t fu = flatOf(u);
                            std::int64_t fv = flatOf(v);
                            if (fu > fv) std::swap(fu, fv);
                            edgePairs.emplace_back(fu, fv);
                        }
                    }
                }
            }
        }
    }
    std::sort(edgePairs.begin(), edgePairs.end());
    edgePairs.erase(std::unique(edgePairs.begin(), edgePairs.end()), edgePairs.end());

    // Points: ascending lattice order; coordinates use the SAME
    // expression as whichever serial mesh rank 0 assembles (see
    // Refine.hpp) so gathered values are bit-identical to the serial
    // classification's inputs.
    const bool baseFormula = regions.empty();
    auto latticeLess = [](const Key& u, const Key& v) {
        if (u.k != v.k) return u.k < v.k;
        if (u.j != v.j) return u.j < v.j;
        return u.i < v.i;
    };
    std::sort(pointKeys.begin(), pointKeys.end(), latticeLess);
    pointKeys.erase(std::unique(pointKeys.begin(), pointKeys.end()), pointKeys.end());
    for (const Key& pk : pointKeys) {
        const int i = pk.i;
        const int j = pk.j;
        const int k = pk.k;
        Vec3 pt;
        if (baseFormula) {
            pt = Vec3{cfg.min.x + (cfg.max.x - cfg.min.x) * (static_cast<double>(i) / nx),
                      cfg.min.y + (cfg.max.y - cfg.min.y) * (static_cast<double>(j) / ny),
                      cfg.min.z + (cfg.max.z - cfg.min.z) * (static_cast<double>(k) / nz)};
        } else {
            pt = Vec3{cfg.min.x + (cfg.max.x - cfg.min.x) * static_cast<double>(i) / (nx * R),
                      cfg.min.y + (cfg.max.y - cfg.min.y) * static_cast<double>(j) / (ny * R),
                      cfg.min.z + (cfg.max.z - cfg.min.z) * static_cast<double>(k) / (nz * R)};
        }
        out.points.push_back(pt);
        out.pointFlat.push_back(fineLatticeFlat(i, j, k, nxF, nyF));
        // The 8 fine cells around the lattice point; a probe outside the
        // domain reads -1 and is skipped. Owned-face points sit at most one
        // fine cell from the owned box, inside the exchanged halo.
        int finest = 0;
        for (int dk = -1; dk <= 0; ++dk) {
            for (int dj = -1; dj <= 0; ++dj) {
                for (int di = -1; di <= 0; ++di) {
                    const int id = tree.at(i + di, j + dj, k + dk);
                    if (id >= 0) finest = std::max(finest, tree.leaf(id).level);
                }
            }
        }
        out.pointFinestLevel.push_back(finest);
    }
    auto localOf = [&](std::int64_t f) {
        return static_cast<int>(std::lower_bound(out.pointFlat.begin(), out.pointFlat.end(), f) -
                                out.pointFlat.begin());
    };
    out.edges.reserve(edgePairs.size());
    for (const auto& [fu, fv] : edgePairs) {
        out.edges.push_back(makeEdgeKey(localOf(fu), localOf(fv)));
    }
    return out;
}

AssembledRefine assembleRefined(const MeshConfig& cfg, const std::vector<LeafRec>& leaves, int R) {
    const int nx = cfg.nx;
    const int ny = cfg.ny;
    const int nz = cfg.nz;
    const int nxF = nx * R;
    const int nyF = ny * R;
    const int nzF = nz * R;
    int maxLevel = 0;
    {
        int r = R;
        while (r > 1) {
            r >>= 1;
            ++maxLevel;
        }
    }

    // Global fine-cell -> leaf lookup, rebuilt from the gathered leaf
    // set (only equality of the returned ids is used below).
    const int baseLo[3] = {0, 0, 0};
    const int baseHi[3] = {nx, ny, nz};
    LocalOctree lookup(baseLo, baseHi, R);
    for (const LeafRec& L : leaves) lookup.refineTo(L.i0, L.j0, L.k0, L.level);
    auto at = [&](int i, int j, int k) { return lookup.at(i, j, k); };

    // --- From here on: the serial assembly, identical in
    // logic (deterministic leaf numbering; quad emission; generalized
    // conformance splice; sorted-key point ids; signature adjacency;
    // upper-triangular internal faces; patches in dict order) -- a pure
    // function of the final leaf set, so byte-identical for every rank
    // count. Only the container changed (the gathered LeafRec list
    // rebuilt into a sparse LocalOctree for the neighbour probes).
    std::vector<int> liveIds(leaves.size());
    for (std::size_t i = 0; i < leaves.size(); ++i) liveIds[i] = static_cast<int>(i);
    std::sort(liveIds.begin(), liveIds.end(), [&](int a, int b) {
        const LeafRec& la = leaves[static_cast<std::size_t>(a)];
        const LeafRec& lb = leaves[static_cast<std::size_t>(b)];
        const int baseA = baseCellIndex(la.i0 / R, la.j0 / R, la.k0 / R, nx, ny);
        const int baseB = baseCellIndex(lb.i0 / R, lb.j0 / R, lb.k0 / R, nx, ny);
        if (baseA != baseB) return baseA < baseB;
        if (la.k0 != lb.k0) return la.k0 < lb.k0;
        if (la.j0 != lb.j0) return la.j0 < lb.j0;
        return la.i0 < lb.i0;
    });
    const int totalCells = static_cast<int>(liveIds.size());
    // One patch per distinct NAME, not per side: several domain sides
    // may share a name and OpenFOAM rejects a repeated patch name in
    // constant/polyMesh/boundary (see mergePatchSpecs). Mirrors
    // BaseMesh's own boundary assembly.
    const MergedPatches merged = mergePatchSpecs(cfg.patches);
    const std::unordered_map<std::string, std::size_t>& sideKeyToPatchOrdinal = merged.sideKeyToOrdinal;

    // Pass A: generate every final leaf's outward quads, in final-id
    // order.
    QuadStore rawQuads;
    rawQuads.cell.reserve(static_cast<std::size_t>(totalCells) * 6);

    auto emitDir = [&](int I0, int J0, int K0, int step, int F, int emittingId, int sideOrd, int dir) {
        for (int p = 0; p < F; ++p) {
            for (int q = 0; q < F; ++q) {
                const std::array<Key, 4> corners = subQuadCorners(dir, I0, J0, K0, step, F, p, q);
                rawQuads.append(corners.data(), 4, emittingId, sideOrd);
            }
        }
    };

    const char* sideKeys6[6] = {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};
    for (int fid = 0; fid < totalCells; ++fid) {
        const LeafRec& L = leaves[static_cast<std::size_t>(liveIds[static_cast<std::size_t>(fid)])];
        const int size = R >> L.level;
        for (int dir = 0; dir < 6; ++dir) {
            int I0 = L.i0, J0 = L.j0, K0 = L.k0;
            int probeI = L.i0, probeJ = L.j0, probeK = L.k0;
            bool isBoundary = false;
            switch (dir) {
                case 0: isBoundary = (L.i0 == 0); probeI = L.i0 - 1; break;
                case 1: isBoundary = (L.i0 + size == nxF); probeI = L.i0 + size; break;
                case 2: isBoundary = (L.j0 == 0); probeJ = L.j0 - 1; break;
                case 3: isBoundary = (L.j0 + size == nyF); probeJ = L.j0 + size; break;
                case 4: isBoundary = (L.k0 == 0); probeK = L.k0 - 1; break;
                case 5: isBoundary = (L.k0 + size == nzF); probeK = L.k0 + size; break;
                default: break;
            }
            if (isBoundary) {
                if (sideKeyToPatchOrdinal.find(sideKeys6[dir]) == sideKeyToPatchOrdinal.end()) {
                    throw std::runtime_error(std::string("Refine error: no patch declared for domain side '") +
                                              sideKeys6[dir] + "'");
                }
                emitDir(I0, J0, K0, size, 1, fid, dir, dir);
                continue;
            }
            // Determine F by sampling the 4 quadrant-origin fine cells on
            // the neighbour side.
            int F = 1;
            if (size >= 2) {
                const int half = size / 2;
                int ids[4];
                int n = 0;
                for (int a = 0; a <= half; a += half) {
                    for (int b = 0; b <= half; b += half) {
                        int qi = probeI, qj = probeJ, qk = probeK;
                        switch (dir) {
                            case 0:
                            case 1:
                                qj = L.j0 + a;
                                qk = L.k0 + b;
                                break;
                            case 2:
                            case 3:
                                qi = L.i0 + a;
                                qk = L.k0 + b;
                                break;
                            default:
                                qi = L.i0 + a;
                                qj = L.j0 + b;
                                break;
                        }
                        ids[n++] = at(qi, qj, qk);
                    }
                }
                const bool allSame = ids[0] == ids[1] && ids[0] == ids[2] && ids[0] == ids[3];
                if (!allSame) F = 2;
            }
            emitDir(I0, J0, K0, size, F, fid, -1, dir);
        }
    }

    // Conformance pass, GENERALIZED: every lattice point strictly inside
    // a quad edge that is used as a corner anywhere is spliced in,
    // sorted along the edge (CornerLines).
    {
        const CornerLines corners(rawQuads.keys);
        QuadStore spliced;
        spliced.cell.reserve(rawQuads.cell.size());
        std::vector<Key> loop;
        const int nQ = rawQuads.size();
        for (int q = 0; q < nQ; ++q) {
            const Key* pts = rawQuads.ptsOf(q);
            const int n = rawQuads.count(q);
            loop.clear();
            for (int e = 0; e < n; ++e) {
                loop.push_back(pts[e]);
                corners.appendInterior(pts[e], pts[(e + 1) % n], loop);
            }
            spliced.append(loop.data(), static_cast<int>(loop.size()), rawQuads.cell[static_cast<std::size_t>(q)],
                           rawQuads.side[static_cast<std::size_t>(q)]);
        }
        rawQuads = std::move(spliced);
    }

    // Collect distinct point Keys, assign ids in sorted order --
    // formerly a std::map<Key,int>; now the sorted-unique vector IS the
    // id assignment (index == id, ascending Key order, identical ids).
    std::vector<Key> sortedKeys(rawQuads.keys);
    std::sort(sortedKeys.begin(), sortedKeys.end());
    sortedKeys.erase(std::unique(sortedKeys.begin(), sortedKeys.end()), sortedKeys.end());
    sortedKeys.shrink_to_fit();
    auto idOfKey = [&](const Key& k) {
        return static_cast<int>(std::lower_bound(sortedKeys.begin(), sortedKeys.end(), k) - sortedKeys.begin());
    };

    // Per-quad point-id loops, flat (same offsets as the spliced quads).
    std::vector<int> quadIdsFlat(rawQuads.keys.size());
    for (std::size_t t = 0; t < rawQuads.keys.size(); ++t) {
        quadIdsFlat[t] = idOfKey(rawQuads.keys[t]);
    }

    // Pass B: adjacency via sorted point-id signature -- formerly an
    // unordered_map keyed by the sorted id vector (one heap node + one
    // heap vector per FACE); now a sort-and-group over quad indices:
    // two quads share a face iff their sorted id signatures are equal,
    // so sorting quad indices by signature makes matching quads
    // adjacent. firstCell/secondCell are assigned in ascending quad
    // order within a group, exactly the order the map registered them.
    std::vector<int> sigFlat(quadIdsFlat);
    {
        const int nQ = rawQuads.size();
        for (int q = 0; q < nQ; ++q) {
            std::sort(sigFlat.begin() + rawQuads.offsets[static_cast<std::size_t>(q)],
                      sigFlat.begin() + rawQuads.offsets[static_cast<std::size_t>(q) + 1]);
        }
    }
    std::vector<int> recFirst(static_cast<std::size_t>(rawQuads.size()), -1);
    std::vector<int> recSecond(static_cast<std::size_t>(rawQuads.size()), -1);
    {
        const int nQ = rawQuads.size();
        std::vector<int> order(static_cast<std::size_t>(nQ));
        for (int q = 0; q < nQ; ++q) order[static_cast<std::size_t>(q)] = q;
        auto sigLess = [&](int a, int b) {
            const int* ab = sigFlat.data() + rawQuads.offsets[static_cast<std::size_t>(a)];
            const int* ae = sigFlat.data() + rawQuads.offsets[static_cast<std::size_t>(a) + 1];
            const int* bb = sigFlat.data() + rawQuads.offsets[static_cast<std::size_t>(b)];
            const int* be = sigFlat.data() + rawQuads.offsets[static_cast<std::size_t>(b) + 1];
            if (std::lexicographical_compare(ab, ae, bb, be)) return true;
            if (std::lexicographical_compare(bb, be, ab, ae)) return false;
            return a < b; // equal signatures: ascending quad order
        };
        std::sort(order.begin(), order.end(), sigLess);
        auto sigEqual = [&](int a, int b) {
            const int* ab = sigFlat.data() + rawQuads.offsets[static_cast<std::size_t>(a)];
            const int* ae = sigFlat.data() + rawQuads.offsets[static_cast<std::size_t>(a) + 1];
            const int* bb = sigFlat.data() + rawQuads.offsets[static_cast<std::size_t>(b)];
            const int* be = sigFlat.data() + rawQuads.offsets[static_cast<std::size_t>(b) + 1];
            return (ae - ab) == (be - bb) && std::equal(ab, ae, bb);
        };
        std::size_t g0 = 0;
        while (g0 < order.size()) {
            std::size_t g1 = g0 + 1;
            while (g1 < order.size() && sigEqual(order[g0], order[g1])) ++g1;
            // Group members are in ascending quad order (the comparator's
            // tie-break); mirror the former map semantics: first
            // registrant fills firstCell, every later one overwrites
            // secondCell.
            int first = rawQuads.cell[static_cast<std::size_t>(order[g0])];
            int second = -1;
            for (std::size_t m = g0 + 1; m < g1; ++m) {
                second = rawQuads.cell[static_cast<std::size_t>(order[m])];
            }
            for (std::size_t m = g0; m < g1; ++m) {
                recFirst[static_cast<std::size_t>(order[m])] = first;
                recSecond[static_cast<std::size_t>(order[m])] = second;
            }
            g0 = g1;
        }
    }
    sigFlat.clear();
    sigFlat.shrink_to_fit();

    // Pass C: emit faces once, from the owner (min cell id) side.
    std::vector<FaceStore> boundaryByPatch(merged.patches.size());
    FaceStore internalFaces;
    {
        const int nQ = rawQuads.size();
        std::vector<int> ids;
        for (int q = 0; q < nQ; ++q) {
            const bool isBoundary = (recSecond[static_cast<std::size_t>(q)] == -1);
            const int owner = isBoundary
                                  ? recFirst[static_cast<std::size_t>(q)]
                                  : std::min(recFirst[static_cast<std::size_t>(q)], recSecond[static_cast<std::size_t>(q)]);
            const int neighbour =
                isBoundary ? -1
                           : std::max(recFirst[static_cast<std::size_t>(q)], recSecond[static_cast<std::size_t>(q)]);
            if (owner != rawQuads.cell[static_cast<std::size_t>(q)]) {
                continue;
            }
            ids.assign(quadIdsFlat.begin() + rawQuads.offsets[static_cast<std::size_t>(q)],
                       quadIdsFlat.begin() + rawQuads.offsets[static_cast<std::size_t>(q) + 1]);
            if (isBoundary) {
                if (rawQuads.side[static_cast<std::size_t>(q)] < 0) {
                    throw std::runtime_error("Refine internal error: unmatched face is not on a domain boundary");
                }
                const std::size_t pOrd =
                    sideKeyToPatchOrdinal.at(sideKeys6[rawQuads.side[static_cast<std::size_t>(q)]]);
                boundaryByPatch[pOrd].append(ids, owner, neighbour, static_cast<int>(pOrd));
            } else {
                internalFaces.append(ids, owner, neighbour, -1);
            }
        }
    }

    internalFaces.stableSortByOwnerNeighbour();

    GeneratedMesh mesh;
    mesh.nInternalFaces = internalFaces.size();
    mesh.faces = std::move(internalFaces);
    for (std::size_t p = 0; p < merged.patches.size(); ++p) {
        PatchInfo info = merged.patches[p];
        info.startFace = mesh.faces.size();
        info.nFaces = boundaryByPatch[p].size();
        mesh.faces.appendAll(boundaryByPatch[p]);
        mesh.patches.push_back(info);
    }

    // Points + pointLevel (+ pointFlat, the key map for the
    // gathered-CutData reconstruction).
    mesh.points.resize(sortedKeys.size());
    std::vector<int> pointLevel(sortedKeys.size(), 0);
    std::vector<long long> pointFlat(sortedKeys.size(), 0);
    auto axisLevel = [&](int v) {
        const int o = v % R;
        if (o == 0) return 0;
        const int tz = trailingZeros(o);
        return maxLevel - tz;
    };
    for (std::size_t id = 0; id < sortedKeys.size(); ++id) {
        const Key& k = sortedKeys[id];
        mesh.points[static_cast<std::size_t>(id)] =
            Vec3{cfg.min.x + (cfg.max.x - cfg.min.x) * static_cast<double>(k.i) / (nx * R),
                 cfg.min.y + (cfg.max.y - cfg.min.y) * static_cast<double>(k.j) / (ny * R),
                 cfg.min.z + (cfg.max.z - cfg.min.z) * static_cast<double>(k.k) / (nz * R)};
        pointLevel[static_cast<std::size_t>(id)] = std::max({axisLevel(k.i), axisLevel(k.j), axisLevel(k.k)});
        pointFlat[static_cast<std::size_t>(id)] = fineLatticeFlat(k.i, k.j, k.k, nxF, nyF);
    }

    // Cells + cellLevel.
    std::vector<int> cellLevel(static_cast<std::size_t>(totalCells), 0);
    for (int fid = 0; fid < totalCells; ++fid) {
        cellLevel[static_cast<std::size_t>(fid)] =
            leaves[static_cast<std::size_t>(liveIds[static_cast<std::size_t>(fid)])].level;
    }
    buildCellFaces(mesh, totalCells);

    AssembledRefine result;
    result.mesh = std::move(mesh);
    result.cellLevel = std::move(cellLevel);
    result.pointLevel = std::move(pointLevel);
    result.pointFlat = std::move(pointFlat);
    return result;
}

RefineStats computeRefineStats(const GeneratedMesh& mesh, const std::vector<int>& cellLevel) {
    RefineStats stats;
    for (int i = 0; i < mesh.nInternalFaces; ++i) {
        const int diff = std::abs(cellLevel[static_cast<std::size_t>(mesh.faces.owner[static_cast<std::size_t>(i)])] -
                                   cellLevel[static_cast<std::size_t>(mesh.faces.neighbour[static_cast<std::size_t>(i)])]);
        stats.maxAdjacentLevelDiff = std::max(stats.maxAdjacentLevelDiff, diff);
    }
    for (int lvl : cellLevel) {
        if (static_cast<std::size_t>(lvl) >= stats.countByLevel.size()) {
            stats.countByLevel.resize(static_cast<std::size_t>(lvl) + 1, 0);
        }
        ++stats.countByLevel[static_cast<std::size_t>(lvl)];
    }
    return stats;
}

namespace {

using CoordKey = std::tuple<double, double, double>;
CoordKey coordKey(const Vec3& v) { return {v.x, v.y, v.z}; }

std::vector<int> uniquePointsOfCell(const GeneratedMesh& mesh, int cellIdx) {
    std::set<int> pts;
    for (int f : mesh.cellFacesOf(cellIdx)) {
        for (int p : mesh.faces.pointsOf(f)) {
            pts.insert(p);
        }
    }
    return std::vector<int>(pts.begin(), pts.end());
}

} // namespace

PostCutLevels propagateLevelsThroughCut(const GeneratedMesh& preCut, const std::vector<int>& preCutCellLevel,
                                         const std::vector<int>& preCutPointLevel, const CutData& cd,
                                         const GeneratedMesh& cut, const std::vector<int>* originCell) {
    PostCutLevels result;

    // --- Point levels: exact-coordinate map from the pre-cut mesh's
    // own points, plus cut-intercept points (max of the cut edge's two
    // endpoint levels).
    std::map<CoordKey, int> coordToLevel;
    for (std::size_t i = 0; i < preCut.points.size(); ++i) {
        coordToLevel[coordKey(preCut.points[i])] = preCutPointLevel[i];
    }
    for (const auto& [key, intercept] : cd.edgeIntercept) {
        const int lvl = std::max(preCutPointLevel[static_cast<std::size_t>(key.first)],
                                  preCutPointLevel[static_cast<std::size_t>(key.second)]);
        coordToLevel[coordKey(intercept.point)] = lvl;
    }
    result.pointLevel.resize(cut.points.size());
    for (std::size_t i = 0; i < cut.points.size(); ++i) {
        auto it = coordToLevel.find(coordKey(cut.points[i]));
        if (it == coordToLevel.end()) {
            throw std::runtime_error("Refine error: post-cut point has no traceable level (unexpected)");
        }
        result.pointLevel[i] = it->second;
    }

    // --- Cell levels.
    result.cellLevel.resize(static_cast<std::size_t>(cut.nCells()));
    if (originCell != nullptr) {
        // Fast path: the cutter already threaded each output
        // cell's originating pre-cut cell index through -- O(1) lookup,
        // no search. The AABB-scan fallback below is O(n_out
        // * n_precut) and measured well past bm_scale_refined's
        // 10-minute gate at 1M-cell scale.
        if (originCell->size() != static_cast<std::size_t>(cut.nCells())) {
            throw std::runtime_error("Refine error: originCell size does not match cut mesh cell count");
        }
        for (std::size_t c = 0; c < static_cast<std::size_t>(cut.nCells()); ++c) {
            const int pre = (*originCell)[c];
            if (pre < 0 || static_cast<std::size_t>(pre) >= preCutCellLevel.size()) {
                throw std::runtime_error("Refine error: originCell entry out of range");
            }
            result.cellLevel[c] = preCutCellLevel[static_cast<std::size_t>(pre)];
        }
        return result;
    }

    // --- Fallback: match each kept output cell to its originating
    // pre-cut cell by exact AABB containment of the output cell's
    // centroid (cutMesh never merges/splits cells, and cut cells stay
    // strictly inside their parent's bounds). O(n_out * n_precut) --
    // retained only for callers that don't have an originCell available.
    struct AABB {
        Vec3 lo, hi;
    };
    std::vector<AABB> preBox(static_cast<std::size_t>(preCut.nCells()));
    for (std::size_t c = 0; c < static_cast<std::size_t>(preCut.nCells()); ++c) {
        std::vector<int> pts = uniquePointsOfCell(preCut, static_cast<int>(c));
        AABB box{Vec3{1e300, 1e300, 1e300}, Vec3{-1e300, -1e300, -1e300}};
        for (int p : pts) {
            const Vec3& v = preCut.points[static_cast<std::size_t>(p)];
            box.lo.x = std::min(box.lo.x, v.x);
            box.lo.y = std::min(box.lo.y, v.y);
            box.lo.z = std::min(box.lo.z, v.z);
            box.hi.x = std::max(box.hi.x, v.x);
            box.hi.y = std::max(box.hi.y, v.y);
            box.hi.z = std::max(box.hi.z, v.z);
        }
        preBox[c] = box;
    }

    for (std::size_t c = 0; c < static_cast<std::size_t>(cut.nCells()); ++c) {
        std::vector<int> pts = uniquePointsOfCell(cut, static_cast<int>(c));
        Vec3 centroid{0, 0, 0};
        for (int p : pts) {
            centroid = centroid + cut.points[static_cast<std::size_t>(p)];
        }
        centroid = centroid * (1.0 / static_cast<double>(pts.size()));

        int best = -1;
        double bestDist = std::numeric_limits<double>::max();
        const double eps = 1e-9;
        for (std::size_t pc = 0; pc < preBox.size(); ++pc) {
            const AABB& box = preBox[pc];
            const bool inside = centroid.x >= box.lo.x - eps && centroid.x <= box.hi.x + eps &&
                                 centroid.y >= box.lo.y - eps && centroid.y <= box.hi.y + eps &&
                                 centroid.z >= box.lo.z - eps && centroid.z <= box.hi.z + eps;
            if (!inside) continue;
            const Vec3 mid = (box.lo + box.hi) * 0.5;
            const Vec3 d = mid - centroid;
            const double dist = dot(d, d);
            if (dist < bestDist) {
                bestDist = dist;
                best = static_cast<int>(pc);
            }
        }
        if (best == -1) {
            throw std::runtime_error("Refine error: post-cut cell centroid matched no pre-cut cell AABB");
        }
        result.cellLevel[c] = preCutCellLevel[static_cast<std::size_t>(best)];
    }

    return result;
}

} // namespace ninja
