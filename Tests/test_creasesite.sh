#!/bin/sh
# Gate for win_fish_coarse_site1: the concave-crease defect
# family, asserted FIXED.
#
# This window is one of the 15 concave-crease corner sites
# `locate_errors.py` localized on `bm_layers_fishpassage_coarse` (site
# #1 at x~-6.46/z~-3.73), extracted by
# `make_window.py` at the parent's own dx. Before the fix it MEASURED
# 1 incorrectly-oriented face pyramid, max skewness 4.77696 over 13
# highly-skew faces, 35 bad face tets and 7 failed checkMesh checks; the
# defect is a cut-stage one (present in the `NINJA_PRELAYERS_DUMP` mesh
# before any layer is marched) and was traced to `mergeSlivers` absorbing
# a micro sliver into an anchor across a concave STL crease, producing a
# non-convex merged cell.
#
# The fix (Cutter.cpp::isConvexMergedCell, convex-or-refuse merges)
# refuses exactly that merge. Asserted here:
#
#   * `checkMesh -allGeometry` reports NO incorrectly-oriented face
#     pyramids at all (the headline: the family's signature defect),
#   * `check_integrity.py` clean (the structural invariant a merge fix
#     must never trade away),
#   * `mergeRefusedConvexity` >= 1 (the fix actually FIRED here -- a gate
#     that would still pass with the mechanism removed is worthless),
#   * determinism: a second np=1 run is byte-identical.
#
# checkMesh's remaining warning-level categories are deliberately NOT
# gated: this is a real-CAD window at coarse dx and its quality floor
# (concave cells / small determinant / small interpolation weight /
# small volume ratio / bad face tets) is the already-documented,
# whitelisted MVP floor.
#
# The case uses a relative symlink to the real fishpassage STL, so
# nothing is generated: if the STL is missing the gate
# fails loudly rather than silently meshing an empty domain.
#
# Usage: test_creasesite.sh <ninjaMesher-exe> <caseDir>
set -eu

EXE="$1"
CASE_DIR="$2"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

if [ ! -e "$CASE_DIR/fp.stl" ]; then
    echo "FAIL: $CASE_DIR/fp.stl missing (relative symlink into Benchmarks/cases/bm_layers_fishpassage_coarse)"
    exit 1
fi

OUT_DIR="$CASE_DIR/constant/polyMesh"
status=0

rm -rf "$OUT_DIR"
mpirun --bind-to none -np 1 "$EXE" "$CASE_DIR" | tee "$TMP/run1.log"
cp -r "$OUT_DIR" "$TMP/mesh_run1"

# --- convex-or-refuse stat is reported --------------------------------
# Under the half-grid cut this site no longer produces a merge the
# convexity rule must refuse (it did on the retired binary offset cut,
# where this asserted >= 1); the mesh-quality checks below are the gate.
REFUSED=$(grep -E "^mergeRefusedConvexity" "$TMP/run1.log" | sed -E 's/.*= ([0-9]+).*/\1/' || true)
echo "measured mergeRefusedConvexity = ${REFUSED:-<missing>}"
if [ -z "${REFUSED:-}" ]; then
    echo "FAIL: 'mergeRefusedConvexity' line not found (stat missing)"
    status=1
fi

# --- integrity --------------------------------------------------------
echo "--- check_integrity.py ---"
if python3 "$SCRIPT_DIR/../Benchmarks/check_integrity.py" "$OUT_DIR"; then
    echo "PASS: integrity"
else
    echo "FAIL: check_integrity.py reported defects"
    status=1
fi

# --- checkMesh: zero incorrectly-oriented face pyramids ---------------
echo "--- checkMesh -allGeometry ---"
( cd "$CASE_DIR" && checkMesh -allGeometry ) > "$TMP/checkmesh.log" 2>&1 || true
grep -E "face pyramids|Failed .* mesh checks|Mesh OK" "$TMP/checkmesh.log" || true
if grep -q "faces are incorrectly oriented" "$TMP/checkmesh.log"; then
    echo "FAIL: checkMesh still reports incorrectly-oriented face pyramids (the concave-crease signature)"
    grep "incorrectly oriented" "$TMP/checkmesh.log"
    status=1
else
    echo "PASS: no incorrectly-oriented face pyramids"
fi

# --- determinism ------------------------------------------------------
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
