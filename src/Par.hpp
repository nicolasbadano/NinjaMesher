// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#pragma once

// MPI multi-rank: the parallel glue, in ONE place.
//
// Model: proper DOMAIN DECOMPOSITION, not replicated topology. The
// base grid is split into a static 3D block per rank (parDims3 /
// parBlockRange); each rank holds only its own block plus a one-base-
// cell halo and builds a rank-local octree, refines it (halo-exchanged
// 2:1 grading, iterated to an Allreduce fixpoint), enumerates its
// OWNED cut entities and classifies/intercepts them locally. Those
// block-local results are then GATHERED to rank 0 as POD records,
// where the mesh is reassembled by the unchanged serial assembly and
// the serial tail runs (cutMesh -> levels -> mergeSlivers -> layers ->
// dropDisconnectedCells -> writer). This is the memory-scaling
// arrangement: no rank except root ever holds the full mesh.
//
// np-invariance rests on the A+parity design's core property that
// every shared quantity is a pure function of (global position,
// replicated STL) -- so any partition of the work produces
// bit-identical values, and assembly is trivially canonical. Boundary
// entities computed by more than one rank are asserted consistent at
// the gather rather than exchanged.
//
// (An earlier attempt used "replicated global topology + partitioned
// work"; it was np-invariant but had NO memory scaling and was
// rejected. Only this collectives module survived from it.)
//
// Every helper here is safe to call when MPI is NOT initialized (plain
// serial test binaries): it then behaves as rank 0 of a 1-rank world,
// and the gathers are no-ops. There is deliberately no `if (np == 1)`
// anywhere else in the codebase -- np=1 runs the exact same code with
// trivial (full-range) partitions and self-gathers.

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

namespace ninja {

// Rank / size of MPI_COMM_WORLD (0 / 1 when MPI is not initialized).
int parRank();
int parSize();

// Contiguous block [lo, hi) of `n` items owned by rank `r` of `p`
// (deterministic; the first n%p ranks get one extra item; empty ranges
// for r >= n are fine -- the degenerate-partition cases).
void parBlockRange(std::size_t n, int r, int p, std::size_t& lo, std::size_t& hi);

// This rank's own block of `n` items.
inline void parMyRange(std::size_t n, std::size_t& lo, std::size_t& hi) {
    parBlockRange(n, parRank(), parSize(), lo, hi);
}

// In-place Allgatherv of a byte array whose per-rank slices follow
// parBlockRange(data.size(), r, p): each rank must have filled its own
// slice; afterwards every rank holds the full array. No-op serially.
void parAllgatherSlices(std::vector<std::uint8_t>& data);

// Allgatherv of raw bytes; returns the concatenation in rank order.
std::vector<std::uint8_t> parAllgatherBytes(const void* local, std::size_t nBytes);

// Gatherv of raw bytes to rank 0 only (non-root ranks get an empty
// vector) -- the memory-scaling counterpart of parAllgatherBytes for
// the rank-0 gather seam.
std::vector<std::uint8_t> parGatherBytesToRoot(const void* local, std::size_t nBytes);

// Symmetric halo exchange primitive: MPI_Sendrecv of raw bytes.
// `dest`/`src` == -1 means no partner on that side (MPI_PROC_NULL);
// recv buffer must be pre-sized to the expected byte count. No-op
// serially (both partners -1 by construction at np=1).
void parSendRecvBytes(int dest, const void* sendBuf, std::size_t sendBytes, int src, void* recvBuf,
                      std::size_t recvBytes, int tag);

// Elementwise MPI_SUM Allreduce over `n` ints (no-op serially).
void parAllreduceSumInts(int* vals, int n);

// Balanced 3D process grid for parSize() ranks (MPI_Dims_create) --
// the static 3D block decomposition of the base grid.
// {1,1,1} when MPI is not initialized.
std::array<int, 3> parDims3();

// Allgatherv of trivially-copyable records, concatenated in rank order
// (rank 0's records first) -- so a loop partitioned by parMyRange over
// a globally-ordered entity list reassembles in exactly the serial
// evaluation order.
template <typename T>
std::vector<T> parAllgatherRecords(const std::vector<T>& local) {
    static_assert(std::is_trivially_copyable<T>::value, "POD records only");
    std::vector<std::uint8_t> bytes =
        parAllgatherBytes(local.empty() ? nullptr : local.data(), local.size() * sizeof(T));
    std::vector<T> out(bytes.size() / sizeof(T));
    if (!out.empty()) {
        std::memcpy(out.data(), bytes.data(), bytes.size());
    }
    return out;
}

// Gatherv of trivially-copyable records to rank 0 (rank order); other
// ranks receive an empty vector.
template <typename T>
std::vector<T> parGatherRecordsToRoot(const std::vector<T>& local) {
    static_assert(std::is_trivially_copyable<T>::value, "POD records only");
    std::vector<std::uint8_t> bytes =
        parGatherBytesToRoot(local.empty() ? nullptr : local.data(), local.size() * sizeof(T));
    std::vector<T> out(bytes.size() / sizeof(T));
    if (!out.empty()) {
        std::memcpy(out.data(), bytes.data(), bytes.size());
    }
    return out;
}

} // namespace ninja
