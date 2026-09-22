#!/bin/sh
# Gate for the layer march's COVERAGE at refinement-level transitions: the
# number of wall faces that end up without layers, on small cases where
# nothing in the geometry justifies losing them.
#
# Each case localizes one mechanism that dropped layers systematically
# wherever the cell size changed along a wall (all MEASURED on
# bm_layers_wfp, 19474 -> 2161 dropped faces overall):
#   * layers_levelstep_oblique_n6 -- a flat wall 17 degrees off the grid
#     crossing three refinement levels. The half-grid ON band was ONE
#     global scalar sized from the finest level in the case, so on the
#     coarser levels almost every crossing snapped to the edge midpoint (up
#     to h/2 off), the stacks came out at ~0.4 t and the collapsed-prism
#     guard dropped them. 192 dropped faces before, 0 after (per-vertex
#     band = h/4 of the finest leaf touching the vertex).
#   * win_wfp_floor / win_wfp_trans -- grid-aligned sub-boxes of
#     bm_layers_wfp: a coarse floor far from any fine box (270 -> 9) and a
#     level transition on real CAD (685 -> 125). The second also exercises
#     the seam cascade: terrace seam faces were re-extruded, always failed,
#     and one-ring dilated onto the healthy stack beside them every step.
#   * win_gate_seam -- a sub-box of bm_layers_gate (729 -> 374, of which
#     355 were wall faces of a 669-cell fluid component the offset cut
#     disconnects and the flood fill discards: they never march, never ship,
#     and are no longer counted -- 23 on the fluid that is written). With the
#     seams held, a face dropped early could still re-extrude at the
#     landing step beside the stack that had kept marching; the shared
#     front point then marched twice and landed twice on the same wall
#     spot: 10 duplicate-position points and a zero-width crack of
#     back-to-back wall faces. Asserted clean by check_integrity.py (a
#     front point is extruded at most once, Layers.cpp `faceHeld`).
#
# Asserted: droppedFaces <= the bound given, check_integrity.py clean, and
# checkMesh reports no severely non-orthogonal (> 70 degrees) faces. The last
# is the layer gate's non-orthogonality predicate (Layers.cpp
# `faceNonOrthDegOf`): before it win_wfp_trans carried 17 such faces, worst
# 82.8, on prisms grown from a chain of 140 x 6 mm needle wall faces; after,
# worst 68.6 for 4 more dropped faces.
#
# With a 4th argument "disclosure", the run must also print the layer stage's
# offset-cut disconnection disclosure (offsetCutDisconnectedComponents): the
# fluid the offset cut seals off is reported, not discovered later as a
# volume diff. win_gate_seam has one such component (669 cells).
#
# Usage: test_layerdrops.sh <ninjaMesher-exe> <caseDir> <maxDroppedFaces> [disclosure]
set -eu

EXE="$1"
CASE_DIR="$2"
MAX_DROPPED="$3"
WANT_DISCLOSURE="${4:-}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

OUT_DIR="$CASE_DIR/constant/polyMesh"
rm -rf "$OUT_DIR"

LOG="$(mktemp)"
trap 'rm -f "$LOG"' EXIT
mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" "$CASE_DIR" > "$LOG" 2>&1 || {
    cat "$LOG"
    echo "FAIL: ninjaMesher exited non-zero"
    exit 1
}

status=0

DROPPED=$(grep -E "^droppedFaces" "$LOG" | tail -1 | sed -E 's/.*= ([0-9]+).*/\1/' || true)
echo "measured droppedFaces = ${DROPPED:-<missing>} (bound $MAX_DROPPED)"
if [ -z "${DROPPED:-}" ] || [ "$DROPPED" -gt "$MAX_DROPPED" ]; then
    echo "FAIL: droppedFaces = ${DROPPED:-<missing>} exceeds $MAX_DROPPED"
    status=1
fi

if [ "$WANT_DISCLOSURE" = "disclosure" ] && ! grep -q "^offsetCutDisconnectedComponents = " "$LOG"; then
    echo "FAIL: the offset-cut disconnection disclosure is missing"
    status=1
fi

if ! python3 "$REPO_ROOT/Benchmarks/check_integrity.py" "$OUT_DIR"; then
    echo "FAIL: check_integrity.py reported a structural mesh pathology"
    status=1
fi

CHECK_LOG="$(mktemp)"
trap 'rm -f "$LOG" "$CHECK_LOG"' EXIT
( cd "$CASE_DIR" && checkMesh -allGeometry ) > "$CHECK_LOG" 2>&1 || true
grep -E "non-orthogonality Max|severely non-orthogonal" "$CHECK_LOG" || true
if ! grep -q "non-orthogonality Max" "$CHECK_LOG"; then
    echo "FAIL: checkMesh reported no non-orthogonality"
    status=1
elif grep -q "severely non-orthogonal" "$CHECK_LOG"; then
    echo "FAIL: checkMesh reports severely non-orthogonal faces"
    status=1
fi

exit $status
