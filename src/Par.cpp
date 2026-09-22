// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#include "Par.hpp"

#include <mpi.h>

#include <stdexcept>

namespace ninja {

namespace {

bool mpiLive() {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (initialized == 0) {
        return false;
    }
    int finalized = 0;
    MPI_Finalized(&finalized);
    return finalized == 0;
}

} // namespace

int parRank() {
    if (!mpiLive()) return 0;
    int r = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    return r;
}

int parSize() {
    if (!mpiLive()) return 1;
    int s = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &s);
    return s;
}

void parBlockRange(std::size_t n, int r, int p, std::size_t& lo, std::size_t& hi) {
    if (p <= 0) {
        throw std::runtime_error("Par error: nonpositive rank count");
    }
    const std::size_t P = static_cast<std::size_t>(p);
    const std::size_t R = static_cast<std::size_t>(r);
    const std::size_t base = n / P;
    const std::size_t rem = n % P;
    lo = R * base + (R < rem ? R : rem);
    hi = lo + base + (R < rem ? 1 : 0);
}

void parAllgatherSlices(std::vector<std::uint8_t>& data) {
    if (!mpiLive()) return;
    const int p = parSize();
    if (p == 1) {
        // Still the one shared code path -- a 1-rank Allgatherv is a
        // self-copy; skipping the MPI call is a pure micro-optimization
        // with identical semantics (the slice already IS the array).
        return;
    }
    std::vector<int> counts(static_cast<std::size_t>(p));
    std::vector<int> displs(static_cast<std::size_t>(p));
    for (int r = 0; r < p; ++r) {
        std::size_t lo = 0, hi = 0;
        parBlockRange(data.size(), r, p, lo, hi);
        counts[static_cast<std::size_t>(r)] = static_cast<int>(hi - lo);
        displs[static_cast<std::size_t>(r)] = static_cast<int>(lo);
    }
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, data.data(), counts.data(), displs.data(),
                   MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);
}

std::vector<std::uint8_t> parAllgatherBytes(const void* local, std::size_t nBytes) {
    if (!mpiLive() || parSize() == 1) {
        std::vector<std::uint8_t> out(nBytes);
        if (nBytes != 0) {
            std::memcpy(out.data(), local, nBytes);
        }
        return out;
    }
    const int p = parSize();
    const int myCount = static_cast<int>(nBytes);
    std::vector<int> counts(static_cast<std::size_t>(p), 0);
    MPI_Allgather(&myCount, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    std::vector<int> displs(static_cast<std::size_t>(p), 0);
    std::size_t total = 0;
    for (int r = 0; r < p; ++r) {
        displs[static_cast<std::size_t>(r)] = static_cast<int>(total);
        total += static_cast<std::size_t>(counts[static_cast<std::size_t>(r)]);
    }
    std::vector<std::uint8_t> out(total);
    MPI_Allgatherv(local, myCount, MPI_UNSIGNED_CHAR, out.data(), counts.data(), displs.data(),
                   MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);
    return out;
}

std::vector<std::uint8_t> parGatherBytesToRoot(const void* local, std::size_t nBytes) {
    if (!mpiLive() || parSize() == 1) {
        std::vector<std::uint8_t> out(nBytes);
        if (nBytes != 0) {
            std::memcpy(out.data(), local, nBytes);
        }
        return out;
    }
    const int p = parSize();
    const int r = parRank();
    const int myCount = static_cast<int>(nBytes);
    std::vector<int> counts(r == 0 ? static_cast<std::size_t>(p) : 0, 0);
    MPI_Gather(&myCount, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<int> displs;
    std::size_t total = 0;
    if (r == 0) {
        displs.assign(static_cast<std::size_t>(p), 0);
        for (int i = 0; i < p; ++i) {
            displs[static_cast<std::size_t>(i)] = static_cast<int>(total);
            total += static_cast<std::size_t>(counts[static_cast<std::size_t>(i)]);
        }
    }
    std::vector<std::uint8_t> out(total);
    MPI_Gatherv(local, myCount, MPI_UNSIGNED_CHAR, out.data(), counts.data(), displs.data(),
                MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    return out;
}

void parSendRecvBytes(int dest, const void* sendBuf, std::size_t sendBytes, int src, void* recvBuf,
                      std::size_t recvBytes, int tag) {
    if (!mpiLive() || parSize() == 1) {
        return; // np=1: no partners exist by construction
    }
    MPI_Sendrecv(sendBuf, static_cast<int>(sendBytes), MPI_UNSIGNED_CHAR,
                 dest < 0 ? MPI_PROC_NULL : dest, tag, recvBuf, static_cast<int>(recvBytes),
                 MPI_UNSIGNED_CHAR, src < 0 ? MPI_PROC_NULL : src, tag, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
}

void parAllreduceSumInts(int* vals, int n) {
    if (!mpiLive() || parSize() == 1) return;
    MPI_Allreduce(MPI_IN_PLACE, vals, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

std::array<int, 3> parDims3() {
    std::array<int, 3> dims{1, 1, 1};
    if (!mpiLive()) return dims;
    int d[3] = {0, 0, 0};
    MPI_Dims_create(parSize(), 3, d);
    dims = {d[0], d[1], d[2]};
    return dims;
}

} // namespace ninja
