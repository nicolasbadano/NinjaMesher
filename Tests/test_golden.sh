#!/bin/sh
# Byte-compares (whitespace-normalized) the generated hex_2x1x1 polyMesh
# against the hand-derived golden files.
# Usage: test_golden.sh <ninjaMesher-exe> <repoRoot>
set -eu

EXE="$1"
ROOT="$2"
CASE_DIR="$ROOT/Tests/cases/hex_2x1x1"
GOLDEN_DIR="$ROOT/Tests/golden/hex_2x1x1"
OUT_DIR="$CASE_DIR/constant/polyMesh"

rm -rf "$OUT_DIR"
mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" "$CASE_DIR"

status=0
for f in points faces owner neighbour boundary cellLevel pointLevel; do
    if [ ! -f "$GOLDEN_DIR/$f" ]; then
        echo "FAIL: missing golden file $GOLDEN_DIR/$f"
        status=1
        continue
    fi
    if ! diff -B -w "$GOLDEN_DIR/$f" "$OUT_DIR/$f" >/tmp/ninjamesher_golden_diff.$$ 2>&1; then
        echo "FAIL: $f differs from golden"
        cat /tmp/ninjamesher_golden_diff.$$
        rm -f /tmp/ninjamesher_golden_diff.$$
        status=1
    else
        rm -f /tmp/ninjamesher_golden_diff.$$
    fi
done

exit $status
