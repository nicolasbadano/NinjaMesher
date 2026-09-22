#!/bin/sh
# Gate for classification-ray robustness on an open
# (non-watertight) shell. `Tests/cases/win_synth_openrim` is a minimal,
# fully-understood synthetic reproducer of a ray-classification bug: a
# tiny open-shell
# box (openbox.stl, -x face deliberately omitted, generated fresh here
# by make_stls.py's `open_box_at` mode) sits far from every grid vertex
# but exactly athwart the first fallback ray direction (+x) from the
# domain's far corner -- pre-fix, that single non-degenerate but
# spurious odd-parity crossing floods the corner vertex as SOLID with
# no real wall nearby, producing a "cut edge reports no STL crossing"
# fatal. Post-fix (multi-ray parity voting in src/CutData.cpp), the
# other 11 fallback directions correctly see no crossing, and majority
# vote (11 fluid vs 1 solid) restores the correct answer.
#
# Asserts:
#   - ninjaMesher exits 0 (must not crash)
#   - --cut-stats reports solidVertices == 0 and cutEdges == 0 (the
#     whole domain is genuinely fluid -- the tiny far box must not
#     flood anything, exactly the measured pre-fix defect)
#
# Usage: test_synthopenrim.sh <ninjaMesher-exe> <caseDir>
set -eu

EXE="$1"; CASE_DIR="$2"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

python3 "$SCRIPT_DIR/scripts/make_stls.py" open_box_at 0.3 -0.08 -0.11 0.5 0.02 0.03 \
    "$CASE_DIR/openbox.stl"

OUT_DIR="$CASE_DIR/constant/polyMesh"
rm -rf "$OUT_DIR"
RUN_LOG=$(mktemp)
trap 'rm -f "$RUN_LOG"' EXIT
mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" "$CASE_DIR" | tee "$RUN_LOG"

status=0

if [ ! -d "$OUT_DIR" ]; then
    echo "FAIL: ninjaMesher did not produce $OUT_DIR"
    exit 1
fi

STATS=$("$EXE" --cut-stats "$CASE_DIR")
echo "$STATS"

SOLID=$(echo "$STATS" | grep -E "^\s*solidVertices" | sed -E 's/.*= ([0-9]+).*/\1/')
CUT=$(echo "$STATS" | grep -E "^\s*cutEdges" | sed -E 's/.*= ([0-9]+).*/\1/')

if [ "$SOLID" != "0" ]; then
    echo "FAIL: expected solidVertices == 0 (the far corner must classify fluid), got $SOLID"
    status=1
fi
if [ "$CUT" != "0" ]; then
    echo "FAIL: expected cutEdges == 0 (no spurious wall from the flood), got $CUT"
    status=1
fi

exit $status
