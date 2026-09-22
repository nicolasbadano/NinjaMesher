#!/usr/bin/env python3
"""Measures the layer-grading geometry of an inflation-layer polyMesh.

Backs the offset_graded gate. Walks wall-normal prism stacks outward
from a boundary patch's faces and prints, per sampled stack, the height
of each successive layer (centroid-to-centroid distance between the
face it entered the cell through and the opposite face), then the
aggregate statistics the gate asserts on:

    stacks <n>
    meanHeights <h0> <h1> ...        (mean over stacks, wall layer first)
    meanRatio <r>                    (mean of h_{j+1}/h_j over stacks/layers)
    monotone <0|1>                   (1 iff every stack is thin-at-wall)

Stdlib only; same Newell-normal/2 area convention as patch_areas.py.
A prism cell's "opposite" face is the unique face of that cell sharing
no point with the face we entered through -- true for any prism whose
top/bottom loops are disjoint, which is exactly what Layers.cpp emits.

Usage: layer_grading.py <polyMeshDir> <patchName> <nLayers> [maxStacks]
Exit 2 on a structural surprise (patch missing, stack shorter than
nLayers, no unique opposite face).
"""
import math
import re
import sys


def parse_list_header(text):
    m = re.search(r"\n(\d+)\n\(\n", text)
    if not m:
        sys.exit("malformed polyMesh list file")
    return int(m.group(1)), text[m.end():]


def read_scalar_list(path):
    with open(path) as f:
        n, body = parse_list_header(f.read())
    lines = body.split("\n")
    return [int(lines[i]) for i in range(n)]


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    mesh_dir = sys.argv[1]
    patch = sys.argv[2]
    n_layers = int(sys.argv[3])
    max_stacks = int(sys.argv[4]) if len(sys.argv) > 4 else 40

    with open(mesh_dir + "/points") as f:
        n_pts, body = parse_list_header(f.read())
    plines = body.split("\n")
    pts = [tuple(float(x) for x in plines[i].strip("()").split()) for i in range(n_pts)]

    with open(mesh_dir + "/faces") as f:
        n_faces, fbody = parse_list_header(f.read())
    flines = fbody.split("\n")
    faces = []
    for i in range(n_faces):
        line = flines[i]
        faces.append([int(x) for x in line[line.index("(") + 1:line.index(")")].split()])

    owner = read_scalar_list(mesh_dir + "/owner")
    neigh = read_scalar_list(mesh_dir + "/neighbour")

    # cell -> face list
    n_cells = max(owner) + 1
    cell_faces = [[] for _ in range(n_cells)]
    for fi, c in enumerate(owner):
        cell_faces[c].append(fi)
    for fi, c in enumerate(neigh):
        cell_faces[c].append(fi)

    with open(mesh_dir + "/boundary") as f:
        bnd = f.read()
    start = count = None
    for m in re.finditer(r"(\w+)\n\s*\{([^}]*)\}", bnd):
        if m.group(1) != patch:
            continue
        sf = re.search(r"startFace\s+(\d+)", m.group(2))
        nf = re.search(r"nFaces\s+(\d+)", m.group(2))
        start, count = int(sf.group(1)), int(nf.group(1))
    if start is None:
        sys.exit(2)

    def centroid(fi):
        idxs = faces[fi]
        n = len(idxs)
        return tuple(sum(pts[p][k] for p in idxs) / n for k in range(3))

    def dist(a, b):
        return math.sqrt(sum((a[k] - b[k]) ** 2 for k in range(3)))

    stride = max(1, count // max_stacks)
    heights = []  # per stack: [h_wall, h_1, ...]
    for wf in range(start, start + count, stride):
        stack = []
        cur_face = wf
        cur_cell = owner[wf]
        ok = True
        for _ in range(n_layers):
            cur_set = set(faces[cur_face])
            opposite = [g for g in cell_faces[cur_cell]
                        if g != cur_face and not (cur_set & set(faces[g]))]
            if len(opposite) != 1:
                ok = False
                break
            top = opposite[0]
            stack.append(dist(centroid(cur_face), centroid(top)))
            # step outward across `top`
            nxt = neigh[top] if (top < len(neigh) and owner[top] == cur_cell) else owner[top]
            if nxt == cur_cell or nxt < 0:
                ok = False
                break
            cur_face, cur_cell = top, nxt
        if ok and len(stack) == n_layers:
            heights.append(stack)

    if not heights:
        sys.exit(2)

    n = len(heights)
    mean_h = [sum(s[j] for s in heights) / n for j in range(n_layers)]
    ratios = [s[j + 1] / s[j] for s in heights for j in range(n_layers - 1) if s[j] > 0]
    monotone = all(all(s[j + 1] > s[j] for j in range(n_layers - 1)) for s in heights)

    print("stacks %d" % n)
    print("meanHeights " + " ".join("%.6g" % h for h in mean_h))
    print("meanRatio %.6g" % (sum(ratios) / len(ratios) if ratios else 0.0))
    print("monotone %d" % (1 if monotone else 0))


if __name__ == "__main__":
    main()
