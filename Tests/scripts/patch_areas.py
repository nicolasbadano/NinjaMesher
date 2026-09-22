#!/usr/bin/env python3
"""Prints '<patchName> <area>' per boundary patch of an OpenFOAM polyMesh.

Per-patch wall-area
attribution check for multi-solid layer cases. Stdlib only, same
Newell-normal/2 area convention the mesher itself uses.

Usage: patch_areas.py <polyMeshDir> [patchName ...]
With patch names given, only those are printed (order preserved);
a requested patch missing from the boundary file is an error (exit 2).
"""
import math
import re
import sys


def parse_list_header(text):
    m = re.search(r"\n(\d+)\n\(\n", text)
    if not m:
        sys.exit("malformed polyMesh list file")
    return int(m.group(1)), text[m.end():]


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    mesh_dir = sys.argv[1]
    wanted = sys.argv[2:]

    with open(mesh_dir + "/points") as f:
        n_pts, body = parse_list_header(f.read())
    lines = body.split("\n")
    pts = [tuple(float(x) for x in lines[i].strip("()").split()) for i in range(n_pts)]

    with open(mesh_dir + "/faces") as f:
        n_faces, fbody = parse_list_header(f.read())
    flines = fbody.split("\n")

    def face_area(fid):
        line = flines[fid]
        idxs = [int(x) for x in line[line.index("(") + 1:line.index(")")].split()]
        nx = ny = nz = 0.0
        n = len(idxs)
        for i in range(n):
            x1, y1, z1 = pts[idxs[i]]
            x2, y2, z2 = pts[idxs[(i + 1) % n]]
            nx += (y1 - y2) * (z1 + z2)
            ny += (z1 - z2) * (x1 + x2)
            nz += (x1 - x2) * (y1 + y2)
        return 0.5 * math.sqrt(nx * nx + ny * ny + nz * nz)

    with open(mesh_dir + "/boundary") as f:
        bnd = f.read()
    areas = {}
    order = []
    for m in re.finditer(r"(\w+)\n\s*\{([^}]*)\}", bnd):
        name, block = m.group(1), m.group(2)
        sf = re.search(r"startFace\s+(\d+)", block)
        nf = re.search(r"nFaces\s+(\d+)", block)
        if not sf or not nf:
            continue
        start, count = int(sf.group(1)), int(nf.group(1))
        areas[name] = sum(face_area(i) for i in range(start, start + count))
        order.append(name)

    names = wanted if wanted else order
    for name in names:
        if name not in areas:
            sys.exit(2)
        print(f"{name} {areas[name]:.6f}")


if __name__ == "__main__":
    main()
