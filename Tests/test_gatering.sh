#!/bin/sh
# Regression gate for the "gate-ring" defect family of bm_layers_gate:
# a ring of INVERTED cells along the downstream skinplate of the radial
# gates, at the exactly-grid-aligned plane z = 1.625.
#
# ATTRIBUTION (measured, this window): the family is produced by
# mergeSlivers, NOT by the layer march -- it is present, at identical
# positions and counts, in the NINJA_PRELAYERS_DUMP mesh before a single
# prism is extruded. Merging a sub-threshold sliver into a kept
# neighbour moves the merged cell's centre; on a 0.198 m plate this mesh
# resolves with 1.6 cells the merged polyhedron can end up with its
# centre on the wrong side of one of its own faces, which is exactly
# what checkMesh reports as "face pyramids incorrectly oriented" (and,
# for the same faces, twice over as "Error in face tets").
#
# The fix is the same post-condition the cut applies to its own cells:
# a merge whose merged polyhedron fails checkMesh's per-face pyramid /
# tet-decomposition tests (Cutter.cpp::cellGeometryIsValid) is REFUSED,
# and the sliver stays as its own already-validated cell.
#
# MEASURED on this window (grid-aligned sub-box of bm_layers_gate around
# one gate section, ~6 s), on the PRE-LAYERS mesh:
#     before: 46 faces incorrectly oriented, 92 bad face tets,
#             max skewness 4.30445, 1187 merged slivers, 0 refused
#     after :  0,  0, max skewness 3.71797, 1141 merged, 46 refused
#              on geometry; 22150 cells written either way
#
# Asserted, on the pre-layers (core) mesh AND on the final mesh:
#   * no misoriented face pyramids, no bad face tets, no highly skew
#     faces (the remaining categories are the documented whitelisted
#     quality floor of a coarse real-CAD cut and are not gated),
#   * check_integrity.py clean,
#   * mergeRefusedGeometry >= 1 (the refusal really fires here),
#   * a kept-cell floor, so "remove everything" cannot pass.
#
# Usage: test_gatering.sh <ninjaMesher-exe> <caseDir>
set -eu

EXE="$1"
CASE_DIR="$2"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

OUT_DIR="$CASE_DIR/constant/polyMesh"
rm -rf "$OUT_DIR"
status=0

NINJA_CORE_INTEGRITY=1 NINJA_PRELAYERS_DUMP="$TMP/pre" \
    mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" "$CASE_DIR" 2>&1 | tee "$TMP/run.log"

# The in-pipeline edge-closure contract, at every stage that reports it.
if grep '^coreIntegrity ' "$TMP/run.log" | grep -v 'nonClosedCells=0 *$' > /dev/null; then
    echo "FAIL: a pipeline stage emitted non-edge-closed cells:"
    grep '^coreIntegrity ' "$TMP/run.log" | grep -v 'nonClosedCells=0 *$'
    status=1
fi

MRG=$(sed -n 's/^mergeRefusedGeometry = \([0-9]*\).*/\1/p' "$TMP/run.log" | tail -1)
echo "measured mergeRefusedGeometry = ${MRG:-<missing>}"
if [ -z "${MRG:-}" ] || [ "$MRG" -lt 1 ]; then
    echo "FAIL: the merge geometry post-condition did not fire (mergeRefusedGeometry = ${MRG:-<missing>})"
    status=1
fi

NCELLS=$(sed -n 's/.*wrote [0-9]* points, \([0-9]*\) cells.*/\1/p' "$TMP/run.log" | tail -1)
echo "keptCells=$NCELLS"
if [ -z "$NCELLS" ] || [ "$NCELLS" -lt 21000 ]; then
    echo "FAIL: only $NCELLS cells kept (>= 21000 expected; 22150 measured at the fix)"
    status=1
fi

echo "--- check_integrity.py ---"
if ! python3 "$REPO_ROOT/Benchmarks/check_integrity.py" "$OUT_DIR"; then
    echo "FAIL: check_integrity.py reported a structural mesh pathology"
    status=1
fi

# checkMesh on BOTH the pre-layers core mesh (where the family lives)
# and the final mesh.
cp -r "$CASE_DIR/system" "$TMP/pre/system"
for stage in pre final; do
    if [ "$stage" = pre ]; then D="$TMP/pre"; else D="$CASE_DIR"; fi
    echo "--- checkMesh -allGeometry ($stage) ---"
    ( cd "$D" && checkMesh -allGeometry ) > "$TMP/cm_$stage.log" 2>&1 || true
    grep -E "face pyramids|face tets|Max skewness" "$TMP/cm_$stage.log" || true
    for pat in "incorrectly oriented" "Error in face tets" "highly skew faces"; do
        if grep -q "$pat" "$TMP/cm_$stage.log"; then
            echo "FAIL: checkMesh ($stage) reports '$pat'"
            status=1
        fi
    done
done

if [ "$status" -eq 0 ]; then
    echo "PASS: no misoriented pyramids / bad face tets / highly skew faces, core or final"
fi
exit $status
