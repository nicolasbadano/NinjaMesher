#!/bin/sh
# Gate for win_wfp_apron: the layer march's POST-CONDITION -- it must
# never emit an open cell, a misoriented face pyramid, a bad face tet or
# a highly skew face; where a prism cannot be built validly its face is
# reverted and the mesh stays away from the wall.
#
# The window is a grid-aligned sub-box of bm_layers_wfp (see the case
# dict header). Before the fix it MEASURED 15 open cells, 20 misoriented
# pyramids, 290 bad face tets, max skewness 1133 and a non-converging
# layer-quality gate. Root cause (Layers.cpp): side quads were re-wound
# GEOMETRICALLY (flip if the Newell normal pointed at the centroid) in
# prismValid, in the quality evaluator and in the emitter, so a folded
# quad passed every guard and was written with a winding inconsistent
# with its neighbours (open cell + misoriented pyramid); and the gate's
# re-march loop shipped whatever was still bad at its pass cap.
#
# Asserted:
#   * checkMesh -allGeometry: no open cells, no misoriented face
#     pyramids, no bad face tets, no highly skew faces, no high aspect
#     ratio cells (the remaining categories are the documented
#     whitelisted floor of a coarse real-CAD cut and are not gated),
#   * check_integrity.py clean,
#   * qualityDroppedFaces >= 1 (the in-march guard actually fired here),
#   * the gate converged (no "stacks still bad" warning),
#   * determinism: a second np=1 run is byte-identical.
#
# Usage: test_wfpapron.sh <ninjaMesher-exe> <caseDir>
set -eu

EXE="$1"
CASE_DIR="$2"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

if [ ! -e "$CASE_DIR/wfp_fw.stl" ] || [ ! -e "$CASE_DIR/wfp_ib.stl" ]; then
    echo "FAIL: $CASE_DIR/wfp_*.stl missing (relative symlinks into Benchmarks/cases/bm_layers_wfp)"
    exit 1
fi

OUT_DIR="$CASE_DIR/constant/polyMesh"
status=0

rm -rf "$OUT_DIR"
mpirun --bind-to none -np 1 "$EXE" "$CASE_DIR" 2>&1 | tee "$TMP/run1.log"
cp -r "$OUT_DIR" "$TMP/mesh_run1"

QD=$(grep -E "^qualityDroppedFaces" "$TMP/run1.log" | sed -E 's/.*= ([0-9]+).*/\1/' || true)
echo "measured qualityDroppedFaces = ${QD:-<missing>}"
if [ -z "${QD:-}" ] || [ "$QD" -lt 1 ]; then
    echo "FAIL: in-march quality guard did not fire (qualityDroppedFaces = ${QD:-<missing>})"
    status=1
fi
# COLLAPSED-PRISM guard (the achieved-vs-requested layer height floor).
# This window reproduces the bm_solver_wfp_interfoam finding exactly:
# BEFORE the guard it emitted wall prisms at minLayerHeightFrac 0.0503 --
# a first layer 5% of the thickness the dict asked for, on a mesh
# checkMesh calls "Cell volumes OK". The floor is kMinAchievedHeightFrac
# = 0.35 in Layers.cpp; asserted a hair under it so the gate binds on the
# guard being ACTIVE, not on a floating-point equality.
MLH=$(grep -E "^minLayerHeightFrac" "$TMP/run1.log" | sed -E 's/.*= ([0-9.eE+-]+).*/\1/' || true)
echo "measured minLayerHeightFrac = ${MLH:-<missing>}"
if [ -z "${MLH:-}" ] || ! python3 -c "import sys; sys.exit(0 if float('${MLH:-0}') >= 0.349 else 1)"; then
    echo "FAIL: collapsed layer prism shipped (minLayerHeightFrac = ${MLH:-<missing>}, floor 0.349)"
    status=1
fi

if grep -q "stacks still bad" "$TMP/run1.log"; then
    echo "FAIL: layer-quality gate did not converge"
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
grep -E "Open cells|face pyramids|face tets|skewness|aspect ratio|Failed .* mesh checks|Mesh OK" "$TMP/checkmesh.log" || true
for pat in "Open cells found" "incorrectly oriented" "Error in face tets" "highly skew faces" "High aspect ratio"; do
    if grep -q "$pat" "$TMP/checkmesh.log"; then
        echo "FAIL: checkMesh reports '$pat'"
        status=1
    fi
done
if [ "$status" -eq 0 ]; then
    echo "PASS: no open cells / misoriented pyramids / bad face tets / skew faces / high aspect cells"
fi

rm -rf "$OUT_DIR"
mpirun --bind-to none -np 1 "$EXE" "$CASE_DIR" >"$TMP/run2.log" 2>&1
cp -r "$OUT_DIR" "$TMP/mesh_run2"
if diff -r "$TMP/mesh_run1" "$TMP/mesh_run2" >/dev/null 2>&1; then
    echo "PASS: determinism (two np=1 runs byte-identical)"
else
    echo "FAIL: two np=1 runs produced DIFFERENT polyMesh output"
    status=1
fi

exit $status
