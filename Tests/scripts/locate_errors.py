#!/usr/bin/env python3
"""Localize checkMesh/integrity errors to grid-aligned windows.

Given a case dir that has already been meshed and run
through `checkMesh -allGeometry` (so constant/polyMesh/sets/* exist),
parse the sets that name offending faces/cells (wrongOrientedFaces,
skewFaces, underdeterminedCells, nonClosedCells, concaveCells, etc --
whatever checkMesh actually wrote for this run) plus
check_integrity.py's non-edge-closed cell ids, and print for each
offending face/cell its centroid and a surrounding bounding box, SNAPPED
OUTWARD to the case's own base-mesh lattice (domain min + n from
ninjaMeshDict, uniform cell size dx = (max-min)/n), padded by N base
cells (default 4).

Stdlib only (Tests/scripts/patch_areas.py conventions: hand-rolled
OpenFOAM list-file parsing, no external deps).

Usage:
  locate_errors.py <caseDir> [--pad N] [--sets name1,name2,...]

Prints one block per requested set: set name, count, and either the
individual bounding box of every listed entity (if small) or, more
usefully, the OVERALL bounding box across all entities in that set
(always printed) -- callers pick whichever `make_window.py` should
consume.
"""
import argparse
import re
import sys
from pathlib import Path


def parse_list_header(text):
    m = re.search(r"\n(\d+)\n\(\n", text)
    if not m:
        raise SystemExit("malformed polyMesh list file")
    return int(m.group(1)), text[m.end():]


def read_points(mesh_dir):
    text = (mesh_dir / "points").read_text()
    n, body = parse_list_header(text)
    lines = body.split("\n")
    pts = []
    for i in range(n):
        pts.append(tuple(float(x) for x in lines[i].strip("()").split()))
    return pts


def read_faces(mesh_dir):
    text = (mesh_dir / "faces").read_text()
    n, body = parse_list_header(text)
    lines = body.split("\n")
    faces = []
    for i in range(n):
        line = lines[i]
        idxs = [int(x) for x in line[line.index("(") + 1:line.index(")")].split()]
        faces.append(idxs)
    return faces


def read_owner_neighbour(mesh_dir):
    def read_ints(path):
        text = path.read_text()
        n, body = parse_list_header(text)
        return [int(x) for x in body.split("\n")[:n]]
    return read_ints(mesh_dir / "owner"), read_ints(mesh_dir / "neighbour")


def read_set(sets_dir, name):
    p = sets_dir / name
    if not p.exists():
        return None
    text = p.read_text()
    # cellSet/faceSet/pointSet: header line "class ... Set;" tells kind
    kind = "cellSet"
    m = re.search(r"class\s+(\w+Set)\s*;", text)
    if m:
        kind = m.group(1)
    # OpenFOAM writes small sets inline ("1(165947)") and larger ones
    # multi-line ("26\n(\n165916\n...\n)"); both match count then a
    # parenthesized, whitespace-separated id list.
    m = re.search(r"\n(\d+)\s*\(\s*(.*?)\s*\)\s*\n*// \*+", text, re.S)
    if not m:
        raise SystemExit(f"malformed polyMesh set file: {name}")
    ids = [int(x) for x in m.group(2).split()]
    return kind, ids


def cell_faces(owner, neighbour, ncells):
    cf = {c: [] for c in range(ncells)}
    for fi, o in enumerate(owner):
        cf[o].append(fi)
    for fi, nb in enumerate(neighbour):
        if nb != -1:
            cf[nb].append(fi)
    return cf


def face_centroid(face, pts):
    xs = [pts[p][0] for p in face]
    ys = [pts[p][1] for p in face]
    zs = [pts[p][2] for p in face]
    n = len(face)
    return (sum(xs) / n, sum(ys) / n, sum(zs) / n)


def cell_centroid(cf_faces, faces, pts):
    idxset = set()
    for fi in cf_faces:
        idxset.update(faces[fi])
    xs = [pts[p][0] for p in idxset]
    ys = [pts[p][1] for p in idxset]
    zs = [pts[p][2] for p in idxset]
    n = len(idxset)
    return (sum(xs) / n, sum(ys) / n, sum(zs) / n)


def read_domain(case_dir):
    """Parses the `domain { min (..); max (..); n (..); }` block of
    ninjaMeshDict (stdlib regex, same convention as reading other
    project dicts in Tests/scripts)."""
    dict_path = case_dir / "system" / "ninjaMeshDict"
    text = dict_path.read_text()
    m = re.search(r"domain\s*\{([^}]*)\}", text, re.S)
    block = m.group(1)
    def vec(name):
        mm = re.search(name + r"\s*\(([^)]*)\)", block)
        return tuple(float(x) for x in mm.group(1).split())
    dmin = vec("min")
    dmax = vec("max")
    n = tuple(int(x) for x in re.search(r"(?<![a-zA-Z])n\s*\(([^)]*)\)", block).group(1).split())
    dx = tuple((dmax[i] - dmin[i]) / n[i] for i in range(3))
    return dmin, dmax, n, dx


def snap_bbox(centroids, dmin, dx, pad):
    """Bounding box of `centroids`, snapped OUTWARD to the base-mesh
    lattice (dmin + k*dx), then padded by `pad` base cells on every
    side."""
    xs = [c[0] for c in centroids]
    ys = [c[1] for c in centroids]
    zs = [c[2] for c in centroids]
    lo = [min(xs), min(ys), min(zs)]
    hi = [max(xs), max(ys), max(zs)]
    import math
    bmin = []
    bmax = []
    for i in range(3):
        klo = math.floor((lo[i] - dmin[i]) / dx[i]) - pad
        khi = math.ceil((hi[i] - dmin[i]) / dx[i]) + pad
        bmin.append(dmin[i] + klo * dx[i])
        bmax.append(dmin[i] + khi * dx[i])
    return tuple(bmin), tuple(bmax)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("caseDir")
    ap.add_argument("--pad", type=int, default=4)
    ap.add_argument("--sets", default=None,
                     help="comma-separated set names; default: all sets present "
                          "plus the synthetic 'integrityNonEdgeClosed' cell list")
    args = ap.parse_args()

    case_dir = Path(args.caseDir)
    mesh_dir = case_dir / "constant" / "polyMesh"
    sets_dir = mesh_dir / "sets"

    pts = read_points(mesh_dir)
    faces = read_faces(mesh_dir)
    owner, neighbour = read_owner_neighbour(mesh_dir)
    ncells = max(owner) + 1

    dmin, dmax, n, dx = read_domain(case_dir)
    print(f"domain lattice: min={dmin} dx={dx} n={n}")

    requested = args.sets.split(",") if args.sets else None
    all_names = requested if requested else sorted(p.name for p in sets_dir.glob("*")) if sets_dir.exists() else []

    cf = None  # lazily built (cell->faces), only needed for cell sets

    for name in all_names:
        parsed = read_set(sets_dir, name)
        if parsed is None:
            continue
        kind, ids = parsed
        if not ids:
            continue
        if kind == "faceSet":
            centroids = [face_centroid(faces[fi], pts) for fi in ids]
        elif kind == "cellSet":
            if cf is None:
                cf = cell_faces(owner, neighbour, ncells)
            centroids = [cell_centroid(cf[ci], faces, pts) for ci in ids]
        elif kind == "pointSet":
            centroids = [pts[pi] for pi in ids]
        else:
            continue
        bmin, bmax = snap_bbox(centroids, dmin, dx, args.pad)
        overall_centroid = tuple(sum(c[i] for c in centroids) / len(centroids) for i in range(3))
        print(f"[{name}] kind={kind} count={len(ids)}")
        print(f"  overall centroid = {overall_centroid}")
        print(f"  window bbox (pad={args.pad}) = min{bmin} max{bmax}")
        print(f"  ids (first 10) = {ids[:10]}")

    # Synthetic: integrity's non-edge-closed cells (from check_integrity.py
    # stdout, if the caller passes it via --integrity-cells is future work;
    # here we just recompute it directly since it's cheap).
    from collections import Counter
    cell_faces_map = cell_faces(owner, neighbour, ncells)
    bad_cells = []
    for c in range(ncells):
        edges = Counter()
        for fi in cell_faces_map[c]:
            f = faces[fi]
            m = len(f)
            for a in range(m):
                e = tuple(sorted((f[a], f[(a + 1) % m])))
                edges[e] += 1
        if any(cnt != 2 for cnt in edges.values()):
            bad_cells.append(c)
    if bad_cells:
        centroids = [cell_centroid(cell_faces_map[c], faces, pts) for c in bad_cells]
        bmin, bmax = snap_bbox(centroids, dmin, dx, args.pad)
        overall_centroid = tuple(sum(c[i] for c in centroids) / len(centroids) for i in range(3))
        print(f"[integrityNonEdgeClosed] kind=cellSet count={len(bad_cells)}")
        print(f"  overall centroid = {overall_centroid}")
        print(f"  window bbox (pad={args.pad}) = min{bmin} max{bmax}")
        print(f"  ids = {bad_cells}")


if __name__ == "__main__":
    main()
