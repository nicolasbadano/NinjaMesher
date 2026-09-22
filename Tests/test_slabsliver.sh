#!/bin/sh
# Gate for cut_slabsliver: the sliver threshold is the unlayered cut's
# aspect-ratio guard (aspect <= 1 / volFrac -- see kSmallVolFrac in
# src/Cutter.cpp), and it has no margin to give away.
#
# slab.stl is a block whose top face lies 0.95 to 1.12 thousandths of a cell
# below the grid plane y = 0.3, tilted so the 100 cells of the row under that
# plane sweep slab volume fractions across the threshold (12 triangles;
# regenerate with a plane y = 0.3 - 0.1*(0.00095 + 0.0001 x + 0.00007 z) over
# [-0.5, 1.5]^2, closed down to y = -0.5).
#
# Asserted:
#   * the sweep really straddles the threshold: some of the 100 row cells are
#     removed as slivers and some are kept (measured 18 / 82);
#   * checkMesh -allGeometry reports no high-aspect-ratio cells, and a max
#     aspect ratio above 900 -- so the case is still probing the edge
#     (measured 996.016). Lowering kSmallVolFrac fails the first half,
#     raising it enough to blunt the test fails the second;
#   * check_integrity.py clean.
#
# Usage: test_slabsliver.sh <ninjaMesher-exe> <caseDir>
set -eu

EXE="$1"
CASE_DIR="$2"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

OUT_DIR="$CASE_DIR/constant/polyMesh"
status=0

rm -rf "$OUT_DIR"
mpirun --bind-to none -np 1 "$EXE" "$CASE_DIR" 2>&1 | tee "$TMP/run1.log"

CUT=$(grep -E "^cut: kept " "$TMP/run1.log" | sed -E 's/.*, cut ([0-9]+), removed.*/\1/' | head -1 || true)
echo "measured cut cells = ${CUT:-<missing>}"
if [ -z "${CUT:-}" ]; then
    echo "FAIL: the mesher printed no cut-cell count"
    status=1
elif [ "$CUT" -le 0 ] || [ "$CUT" -ge 100 ]; then
    echo "FAIL: $CUT of the 100 row cells kept -- the sweep no longer straddles the sliver threshold"
    status=1
fi

echo "--- check_integrity.py ---"
if python3 "$SCRIPT_DIR/../Benchmarks/check_integrity.py" "$OUT_DIR"; then
    echo "PASS: integrity"
else
    echo "FAIL: check_integrity.py reported defects"
    status=1
fi

echo "--- checkMesh -allGeometry ---"
( cd "$CASE_DIR" && checkMesh -allGeometry ) > "$TMP/checkmesh.log" 2>&1 || true
grep -E "aspect ratio|Failed .* mesh checks|Mesh OK" "$TMP/checkmesh.log" || true
for pat in "Open cells found" "incorrectly oriented" "Error in face tets" "highly skew faces" "High aspect ratio"; do
    if grep -q "$pat" "$TMP/checkmesh.log"; then
        echo "FAIL: checkMesh reports '$pat'"
        status=1
    fi
done

AR=$(grep -E "Max aspect ratio" "$TMP/checkmesh.log" \
     | sed -E 's/.*Max aspect ratio(:| =) ([0-9.eE+-]+).*/\2/' | head -1 || true)
echo "measured Max aspect ratio = ${AR:-<missing>}"
if [ -z "${AR:-}" ]; then
    echo "FAIL: checkMesh reported no aspect ratio"
    status=1
elif ! python3 -c "import sys; a=float(sys.argv[1]); sys.exit(0 if 900.0 < a <= 1000.0 else 1)" "$AR"; then
    echo "FAIL: max aspect ratio $AR outside (900, 1000]"
    status=1
fi

exit $status
