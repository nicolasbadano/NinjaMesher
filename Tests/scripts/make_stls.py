#!/usr/bin/env python3
"""Generates the ASCII STL fixtures used by the cut-case tests.

Stdlib only. Usage:

    make_stls.py box    <outfile>   # axis-aligned box (0.5,-1,-1)-(3,2,2)
    make_stls.py sphere <outfile>   # icosphere r=0.3, 3 subdivisions, center (0.5,0.5,0.5)
    make_stls.py rotcube <outfile>  # cube side 0.4, 45deg about z, center (0.5,0.5,0.5)
    make_stls.py sphere_at <cx> <cy> <cz> <r> <subdiv> <outfile>
                                     # icosphere of radius r at (cx,cy,cz), <subdiv>
                                     # subdivisions (multi-STL / offcenter cases)
    make_stls.py box_at <x0> <y0> <z0> <x1> <y1> <z1> <outfile>
                                     # axis-aligned box (bm_thin_plate)
    make_stls.py capped_cylinder <cx> <cy> <z0> <z1> <r> <nseg> <outfile>
                                     # closed (capped) cylinder, axis along z,
                                     # from z0 to z1 (bm_int_pipe)
    make_stls.py elbow_tube <x0> <y0> <z0> <x1> <y1> <z1> <x2> <y2> <z2> <r> <nseg> <outfile>
                                     # single watertight compound tube swept
                                     # along the two-segment path
                                     # P0->P1->P2 with a mitred joint at P1
                                     # (bm_int_elbow)
    make_stls.py sphere_f32 <cx> <cy> <cz> <r> <subdiv> <outfile>
                                     # icosphere with every vertex coordinate
                                     # rounded to float32 (see box_f32) --
                                     # many facets means many independent
                                     # near-grid-plane float32 roundings,
                                     # matching the real fishpassage
                                     # multi-facet mechanism that produces
                                     # tiny edges (mixed cells already cut
                                     # by curvature ALSO need a coincident-
                                     # plane intercept nearby).
    make_stls.py cyl_f32 <cy> <cz> <x0> <x1> <r> <nseg> <outfile>
                                     # capped cylinder, axis along x, from
                                     # x0 to x1, float32-rounded (see
                                     # box_f32/sphere_f32 above) -- flat end
                                     # cap at x0 sits on a float32-rounded
                                     # grid-plane coordinate without the
                                     # box's degenerate axis-aligned side-
                                     # wall/cap edge coincidence.
    make_stls.py box_f32 <x0> <y0> <z0> <x1> <y1> <z1> <outfile>
                                     # axis-aligned box, every vertex coordinate
                                     # rounded to the nearest float32 before being
                                     # written out -- simulates what a BINARY STL
                                     # (float32 vertex storage) reader would hand
                                     # NinjaMesher, e.g. an exact grid-plane
                                     # coordinate like 0.3 stored as
                                     # float32(0.3)=0.300000011920929, ~1.19e-8
                                     # off the true double value (the float32
                                     # seam regression case).

Triangles are written CCW as seen from outside the solid (outward
normal); NinjaMesher's ray-parity classification is winding-agnostic
 so exact winding consistency is not load-bearing,
but consistent outward winding is still emitted for good STL hygiene.
"""

import math
import struct
import sys


def to_f32(x):
    """Round a python float (double) to the nearest float32 value,
    returned as a double -- what a binary-STL reader sees when the file
    was written with float32 vertex storage."""
    return struct.unpack("<f", struct.pack("<f", x))[0]


def write_stl(path, name, triangles):
    with open(path, "w") as f:
        f.write(f"solid {name}\n")
        for (v0, v1, v2) in triangles:
            ux, uy, uz = (v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2])
            wx, wy, wz = (v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2])
            nx = uy * wz - uz * wy
            ny = uz * wx - ux * wz
            nz = ux * wy - uy * wx
            norm = math.sqrt(nx * nx + ny * ny + nz * nz) or 1.0
            nx, ny, nz = nx / norm, ny / norm, nz / norm
            f.write(f"  facet normal {nx:.6e} {ny:.6e} {nz:.6e}\n")
            f.write("    outer loop\n")
            for v in (v0, v1, v2):
                f.write(f"      vertex {v[0]:.10e} {v[1]:.10e} {v[2]:.10e}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write(f"endsolid {name}\n")


def box_triangles(lo, hi):
    x0, y0, z0 = lo
    x1, y1, z1 = hi
    p = {
        "000": (x0, y0, z0), "100": (x1, y0, z0), "010": (x0, y1, z0), "110": (x1, y1, z0),
        "001": (x0, y0, z1), "101": (x1, y0, z1), "011": (x0, y1, z1), "111": (x1, y1, z1),
    }
    quads = [
        ("000", "010", "011", "001"),  # -x
        ("100", "101", "111", "110"),  # +x
        ("000", "001", "101", "100"),  # -y
        ("010", "110", "111", "011"),  # +y
        ("000", "100", "110", "010"),  # -z
        ("001", "011", "111", "101"),  # +z
    ]
    tris = []
    for a, b, c, d in quads:
        tris.append((p[a], p[b], p[c]))
        tris.append((p[a], p[c], p[d]))
    return tris


def open_box_triangles(lo, hi, omit="-x"):
    """Same axis-aligned box as `box_triangles`, but with ONE named face
    OMITTED -- an open shell (synthetic reproducer for
    the "ray-parity on an open shell floods a far classification point"
    bug measured on the real motorbike CAD's tire/wheel geometry -- see
    Tests/cases/win_synth_openrim). `omit` is one of
    -x/+x/-y/+y/-z/+z, naming which face to leave out."""
    x0, y0, z0 = lo
    x1, y1, z1 = hi
    p = {
        "000": (x0, y0, z0), "100": (x1, y0, z0), "010": (x0, y1, z0), "110": (x1, y1, z0),
        "001": (x0, y0, z1), "101": (x1, y0, z1), "011": (x0, y1, z1), "111": (x1, y1, z1),
    }
    named_quads = {
        "-x": ("000", "010", "011", "001"),
        "+x": ("100", "101", "111", "110"),
        "-y": ("000", "001", "101", "100"),
        "+y": ("010", "110", "111", "011"),
        "-z": ("000", "100", "110", "010"),
        "+z": ("001", "011", "111", "101"),
    }
    tris = []
    for name, (a, b, c, d) in named_quads.items():
        if name == omit:
            continue
        tris.append((p[a], p[b], p[c]))
        tris.append((p[a], p[c], p[d]))
    return tris


def icosphere_triangles(radius, center, subdivisions):
    t = (1.0 + math.sqrt(5.0)) / 2.0
    verts = [
        (-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0),
        (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t),
        (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1),
    ]

    def normalize(v):
        n = math.sqrt(sum(c * c for c in v))
        return tuple(c / n for c in v)

    verts = [normalize(v) for v in verts]
    faces = [
        (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
        (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
        (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
        (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
    ]

    midpoint_cache = {}

    def midpoint(i, j, vlist):
        key = (min(i, j), max(i, j))
        if key in midpoint_cache:
            return midpoint_cache[key]
        a, b = vlist[i], vlist[j]
        m = normalize(((a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2))
        vlist.append(m)
        idx = len(vlist) - 1
        midpoint_cache[key] = idx
        return idx

    vlist = list(verts)
    for _ in range(subdivisions):
        midpoint_cache.clear()
        new_faces = []
        for (i0, i1, i2) in faces:
            a = midpoint(i0, i1, vlist)
            b = midpoint(i1, i2, vlist)
            c = midpoint(i2, i0, vlist)
            new_faces.append((i0, a, c))
            new_faces.append((i1, b, a))
            new_faces.append((i2, c, b))
            new_faces.append((a, b, c))
        faces = new_faces

    tris = []
    for (i0, i1, i2) in faces:
        v0 = vlist[i0]
        v1 = vlist[i1]
        v2 = vlist[i2]
        p0 = (center[0] + radius * v0[0], center[1] + radius * v0[1], center[2] + radius * v0[2])
        p1 = (center[0] + radius * v1[0], center[1] + radius * v1[1], center[2] + radius * v1[2])
        p2 = (center[0] + radius * v2[0], center[1] + radius * v2[1], center[2] + radius * v2[2])
        tris.append((p0, p1, p2))
    return tris


def rotcube_triangles(side, center, angle_deg):
    h = side / 2.0
    local = [
        (-h, -h, -h), (h, -h, -h), (h, h, -h), (-h, h, -h),
        (-h, -h, h), (h, -h, h), (h, h, h), (-h, h, h),
    ]
    a = math.radians(angle_deg)
    ca, sa = math.cos(a), math.sin(a)

    def rot(p):
        x, y, z = p
        return (x * ca - y * sa, x * sa + y * ca, z)

    verts = [tuple(c + o for c, o in zip(center, rot(p))) for p in local]
    quads = [
        (0, 3, 2, 1),  # -z
        (4, 5, 6, 7),  # +z
        (0, 1, 5, 4),  # -y
        (3, 7, 6, 2),  # +y
        (0, 4, 7, 3),  # -x
        (1, 2, 6, 5),  # +x
    ]
    tris = []
    for a4, b4, c4, d4 in quads:
        tris.append((verts[a4], verts[b4], verts[c4]))
        tris.append((verts[a4], verts[c4], verts[d4]))
    return tris


def capped_cylinder_triangles(cx, cy, z0, z1, r, nseg):
    """Closed (capped) cylinder, axis along z, side + two triangle-fan
    caps -- a single trivially-watertight compound, no CSG (
    bm_int_pipe)."""
    ring = [
        (cx + r * math.cos(2 * math.pi * i / nseg), cy + r * math.sin(2 * math.pi * i / nseg))
        for i in range(nseg)
    ]
    tris = []
    # side quads
    for i in range(nseg):
        j = (i + 1) % nseg
        a = (ring[i][0], ring[i][1], z0)
        b = (ring[j][0], ring[j][1], z0)
        c = (ring[j][0], ring[j][1], z1)
        d = (ring[i][0], ring[i][1], z1)
        tris.append((a, b, c))
        tris.append((a, c, d))
    # caps (fan from centre, outward-facing winding)
    c0 = (cx, cy, z0)
    c1 = (cx, cy, z1)
    for i in range(nseg):
        j = (i + 1) % nseg
        p0 = (ring[i][0], ring[i][1], z0)
        p1 = (ring[j][0], ring[j][1], z0)
        tris.append((c0, p1, p0))  # bottom cap, outward normal -z
        q0 = (ring[i][0], ring[i][1], z1)
        q1 = (ring[j][0], ring[j][1], z1)
        tris.append((c1, q0, q1))  # top cap, outward normal +z
    return tris


def _normalize3(v):
    n = math.sqrt(sum(c * c for c in v))
    return tuple(c / n for c in v)


def _ring_at(centre, axis, r, nseg):
    """A circle of radius r, centred at `centre`, lying in the plane
    perpendicular to `axis` (unit vector)."""
    axis = _normalize3(axis)
    # any vector not parallel to axis
    ref = (1.0, 0.0, 0.0) if abs(axis[0]) < 0.9 else (0.0, 1.0, 0.0)
    u = _normalize3((
        ref[1] * axis[2] - ref[2] * axis[1],
        ref[2] * axis[0] - ref[0] * axis[2],
        ref[0] * axis[1] - ref[1] * axis[0],
    ))
    v = (
        axis[1] * u[2] - axis[2] * u[1],
        axis[2] * u[0] - axis[0] * u[2],
        axis[0] * u[1] - axis[1] * u[0],
    )
    pts = []
    for i in range(nseg):
        a = 2 * math.pi * i / nseg
        pts.append(tuple(
            centre[k] + r * (math.cos(a) * u[k] + math.sin(a) * v[k]) for k in range(3)
        ))
    return pts


def elbow_tube_triangles(p0, p1, p2, r, nseg):
    """Single watertight compound: a circular tube swept along the
    2-segment path p0->p1->p2, with a mitred cross-section ring at the
    corner p1 (bisector of the two segment directions) and end caps at
    p0 and p2. No CSG -- a direct sweep (bm_int_elbow; a faceted
    torus-elbow segment is acceptable here)."""
    d0 = _normalize3(tuple(p1[k] - p0[k] for k in range(3)))
    d1 = _normalize3(tuple(p2[k] - p1[k] for k in range(3)))
    bisector = _normalize3(tuple(d0[k] + d1[k] for k in range(3)))

    ring0 = _ring_at(p0, d0, r, nseg)
    ringM = _ring_at(p1, bisector, r, nseg)
    ring2 = _ring_at(p2, d1, r, nseg)

    tris = []

    def side(ringA, ringB):
        out = []
        for i in range(nseg):
            j = (i + 1) % nseg
            a, b, c, d = ringA[i], ringA[j], ringB[j], ringB[i]
            out.append((a, b, c))
            out.append((a, c, d))
        return out

    tris += side(ring0, ringM)
    tris += side(ringM, ring2)

    # end caps (fans), outward normal = -path direction at p0, +path
    # direction at p2
    for i in range(nseg):
        j = (i + 1) % nseg
        tris.append((p0, ring0[j], ring0[i]))
        tris.append((p2, ring2[i], ring2[j]))
    return tris


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    shape = sys.argv[1]
    if shape == "sphere_at":
        if len(sys.argv) != 8:
            print(__doc__)
            return 1
        cx, cy, cz, r, subdiv, outfile = (
            float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]),
            float(sys.argv[5]), int(sys.argv[6]), sys.argv[7],
        )
        tris = icosphere_triangles(r, (cx, cy, cz), subdiv)
        write_stl(outfile, "sphere_at", tris)
        return 0
    if shape == "box_at":
        if len(sys.argv) != 9:
            print(__doc__)
            return 1
        x0, y0, z0, x1, y1, z1 = (float(v) for v in sys.argv[2:8])
        outfile = sys.argv[8]
        tris = box_triangles((x0, y0, z0), (x1, y1, z1))
        write_stl(outfile, "box_at", tris)
        return 0
    if shape == "open_box_at":
        if len(sys.argv) != 9:
            print(__doc__)
            return 1
        x0, y0, z0, x1, y1, z1 = (float(v) for v in sys.argv[2:8])
        outfile = sys.argv[8]
        tris = open_box_triangles((x0, y0, z0), (x1, y1, z1))
        write_stl(outfile, "open_box_at", tris)
        return 0
    if shape == "sphere_f32":
        if len(sys.argv) != 8:
            print(__doc__)
            return 1
        cx, cy, cz, r, subdiv, outfile = (
            float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]),
            float(sys.argv[5]), int(sys.argv[6]), sys.argv[7],
        )
        tris = icosphere_triangles(r, (cx, cy, cz), subdiv)
        tris = [tuple(tuple(to_f32(c) for c in v) for v in tri) for tri in tris]
        write_stl(outfile, "sphere_f32", tris)
        return 0
    if shape == "cyl_f32":
        if len(sys.argv) != 9:
            print(__doc__)
            return 1
        cy, cz, x0, x1, r = (float(v) for v in sys.argv[2:7])
        nseg = int(sys.argv[7])
        outfile = sys.argv[8]
        # capped_cylinder_triangles' axis is z; permute (x,y,z)->(z,x,y) so
        # the cylinder axis becomes x. Only the two end caps are flat
        # (axis-aligned) faces -- the curved side wall has no rectangular
        # axis-aligned faces, so (unlike a box) it shares no edge with the
        # x=x0 cap along a whole grid line, avoiding a degenerate
        # edge-on-face coincidence with the domain grid while still
        # giving a genuine finite flat cap coincident with a grid plane.
        raw = capped_cylinder_triangles(cy, cz, x0, x1, r, nseg)
        tris = [tuple((v[2], v[0], v[1]) for v in tri) for tri in raw]
        tris = [tuple(tuple(to_f32(c) for c in v) for v in tri) for tri in tris]
        write_stl(outfile, "cyl_f32", tris)
        return 0
    if shape == "box_f32":
        if len(sys.argv) != 9:
            print(__doc__)
            return 1
        x0, y0, z0, x1, y1, z1 = (to_f32(float(v)) for v in sys.argv[2:8])
        outfile = sys.argv[8]
        tris = box_triangles((x0, y0, z0), (x1, y1, z1))
        write_stl(outfile, "box_f32", tris)
        return 0
    if shape == "capped_cylinder":
        if len(sys.argv) != 9:
            print(__doc__)
            return 1
        cx, cy, z0, z1, r = (float(v) for v in sys.argv[2:7])
        nseg = int(sys.argv[7])
        outfile = sys.argv[8]
        tris = capped_cylinder_triangles(cx, cy, z0, z1, r, nseg)
        write_stl(outfile, "capped_cylinder", tris)
        return 0
    if shape == "elbow_tube":
        if len(sys.argv) != 14:
            print(__doc__)
            return 1
        vals = [float(v) for v in sys.argv[2:13]]
        p0 = tuple(vals[0:3])
        p1 = tuple(vals[3:6])
        p2 = tuple(vals[6:9])
        r, nseg = vals[9], int(vals[10])
        outfile = sys.argv[13]
        tris = elbow_tube_triangles(p0, p1, p2, r, nseg)
        write_stl(outfile, "elbow_tube", tris)
        return 0

    outfile = sys.argv[2]
    if shape == "box":
        tris = box_triangles((0.5, -1.0, -1.0), (3.0, 2.0, 2.0))
        write_stl(outfile, "boxwall", tris)
    elif shape == "sphere":
        tris = icosphere_triangles(0.3, (0.5, 0.5, 0.5), 3)
        write_stl(outfile, "sphere", tris)
    elif shape == "rotcube":
        tris = rotcube_triangles(0.4, (0.5, 0.5, 0.5), 45.0)
        write_stl(outfile, "rotcube", tris)
    else:
        print(f"unknown shape '{shape}'", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
