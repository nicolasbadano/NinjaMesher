#!/usr/bin/env python3
"""Structural mesh-integrity checks that checkMesh + volume/area gates
provably cannot catch (the bm_adv_gridpoint lesson: a mesh can be
geometrically perfect and checkMesh-acceptable while topologically
cracked — duplicate-position points made neighbouring faces reference
different indices for the same physical point).

Checks, all hard failures:
  1. duplicate-position points (bit-exact coordinate comparison — the
     writer's dedup contract is exact-match, so ANY coincidence is a bug)
  2. per-cell closure: within each cell, every undirected edge of its
     face set must appear exactly TWICE (the polyhedron is watertight
     index-wise, not just coordinate-wise)
  3. degenerate faces: < 3 points, repeated point indices, or < 3
     distinct positions
  4. duplicate faces: two faces with the identical point-index set
  5. unreferenced points

Usage: check_integrity.py <polyMesh_dir>
Exit 0 = clean; 1 = pathology found (details on stdout).
"""

import re
import sys
from collections import Counter, defaultdict


def read_list(path):
    with open(path, "r", errors="replace") as f:
        txt = f.read()
    m = re.search(r"\n(\d+)\s*\n\(\n(.*)\n\)", txt, re.S)
    if not m:
        raise SystemExit(f"cannot parse OpenFOAM list file: {path}")
    n = int(m.group(1))
    body = m.group(2).split("\n")
    if len(body) != n:
        raise SystemExit(f"{path}: header says {n} entries, found {len(body)}")
    return body


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    d = sys.argv[1].rstrip("/")

    pts = [tuple(float(x) for x in re.findall(r"[-+0-9.eE]+", line))
           for line in read_list(f"{d}/points")]
    faces = [[int(x) for x in re.findall(r"\d+", line)[1:]]
             for line in read_list(f"{d}/faces")]
    owner = [int(line) for line in read_list(f"{d}/owner")]
    neighbour = [int(line) for line in read_list(f"{d}/neighbour")]

    failures = []

    # 1. duplicate-position points
    dup = {k: v for k, v in Counter(pts).items() if v > 1}
    if dup:
        sample = list(dup.items())[:5]
        failures.append(f"duplicate-position points: {len(dup)} positions, e.g. {sample}")

    # 3. degenerate faces
    for i, f in enumerate(faces):
        if len(f) < 3:
            failures.append(f"face {i}: fewer than 3 points")
            continue
        if len(set(f)) != len(f):
            failures.append(f"face {i}: repeated point index")
        if len({pts[p] for p in f}) < 3:
            failures.append(f"face {i}: fewer than 3 distinct positions")

    # 4. duplicate faces (same index set)
    sig_count = Counter(tuple(sorted(f)) for f in faces)
    ndup_faces = sum(1 for v in sig_count.values() if v > 1)
    if ndup_faces:
        failures.append(f"duplicate faces (same point set): {ndup_faces}")

    # 5. unreferenced points
    used = {p for f in faces for p in f}
    unref = len(pts) - len(used)
    if unref:
        failures.append(f"unreferenced points: {unref}")

    # 2. per-cell closure (edge appears exactly twice within each cell)
    ncells = max(owner) + 1 if owner else 0
    cell_faces = defaultdict(list)
    for i, o in enumerate(owner):
        cell_faces[o].append(i)
    for i, nb in enumerate(neighbour):
        cell_faces[nb].append(i)
    open_cells = []
    for c in range(ncells):
        edges = Counter()
        for fi in cell_faces[c]:
            f = faces[fi]
            for a in range(len(f)):
                e = tuple(sorted((f[a], f[(a + 1) % len(f)])))
                edges[e] += 1
        bad = [e for e, cnt in edges.items() if cnt != 2]
        if bad:
            open_cells.append((c, len(bad)))
    if open_cells:
        failures.append(
            f"cells not edge-closed (watertightness violation): {len(open_cells)}"
            f" cells, e.g. {open_cells[:5]} (cell, bad-edge-count)")

    if failures:
        print(f"INTEGRITY FAIL ({d}):")
        for msg in failures:
            print(f"  - {msg}")
        return 1
    print(f"integrity ok: {len(pts)} points, {len(faces)} faces, {ncells} cells")
    return 0


if __name__ == "__main__":
    sys.exit(main())
