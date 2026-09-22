#!/bin/sh
# Regression gate for the input-precision on-surface band (see
# classifyVertices/inputPrecisionBand in src/CutData.cpp): a binary STL's
# float32-stored vertex coordinates are typically off an exact
# grid-vertex-coincident plane by ~1e-7 * |coordinate| -- far above the
# old 1e-9 on-surface tolerance -- so such planes used to produce cut
# edges whose intercept landed a hair from the coincident grid vertex
# (checkMesh "Edges too small"). float32_seam/plate.stl is a capped
# cylinder whose flat end cap sits at float32(0.3), a grid-vertex plane
# of the case's n=(15 15 15) domain (h=0.06) -- built by
# `make_stls.py cyl_f32`, which rounds every vertex coordinate to the
# nearest float32 the way a binary-STL reader would.
#
# Usage: test_float32seam.sh <ninjaMesher-exe> <caseDir>
set -eu

EXE="$1"
CASE_DIR="$2"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

python3 "$SCRIPT_DIR/scripts/make_stls.py" cyl_f32 0.45 0.45 0.3 0.66 0.24 24 "$CASE_DIR/plate.stl"

OUT_DIR="$CASE_DIR/constant/polyMesh"
rm -rf "$OUT_DIR"
mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" "$CASE_DIR"

status=0

# Edge-length gate: 1e-4 is well below the domain's grid size (h=0.06)
# and far above the ~1.2e-8 float32 noise floor the band is designed to
# absorb -- any surviving edge below this threshold means the band
# failed to snap a near-coincident vertex/intercept.
set +e
EDGE_CHECK=$(python3 - "$OUT_DIR" <<'PYEOF'
import math
import re
import sys

d = sys.argv[1]


def read_list(path):
    with open(path) as f:
        txt = f.read()
    m = re.search(r"\n(\d+)\s*\n\(\n(.*)\n\)", txt, re.S)
    n = int(m.group(1))
    body = m.group(2).split("\n")
    assert len(body) == n, f"{path}: header says {n}, found {len(body)}"
    return body


pts = [tuple(float(x) for x in re.findall(r"[-+0-9.eE]+", line))
       for line in read_list(f"{d}/points")]
faces = [[int(x) for x in re.findall(r"\d+", line)[1:]]
         for line in read_list(f"{d}/faces")]

THRESH = 1e-4
tiny = []
min_edge = float("inf")
for fi, f in enumerate(faces):
    n = len(f)
    for i in range(n):
        a, b = pts[f[i]], pts[f[(i + 1) % n]]
        dist = math.dist(a, b)
        if dist > 0:
            min_edge = min(min_edge, dist)
        if 0 < dist < THRESH:
            tiny.append((fi, dist))

print(f"min_edge={min_edge:.6e} n_tiny={len(tiny)}")
if tiny:
    for fi, dist in tiny[:10]:
        print(f"FAIL: face {fi} has edge length {dist:.6e} < {THRESH:.1e}")
    sys.exit(1)
sys.exit(0)
PYEOF
)
EDGE_STATUS=$?
set -e
echo "$EDGE_CHECK"
if [ "$EDGE_STATUS" -ne 0 ]; then
    echo "FAIL: tiny edges present (input-precision band did not absorb the float32 seam)"
    status=1
fi

if ! python3 "$REPO_ROOT/Benchmarks/check_integrity.py" "$OUT_DIR"; then
    echo "FAIL: check_integrity.py reported a structural mesh pathology"
    status=1
fi

exit $status
