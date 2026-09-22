#!/bin/sh
# Regression gate for the 2:1-interface hanging-node inconsistency at an
# exactly-coplanar grid plane (see Tests/cases/win_gate_hanging's dict
# header for the full mechanism and provenance, and src/Cutter.cpp's
# chord-degeneracy branch for the fix).
#
# The contract asserted here is the per-cell edge-closure one: every
# emitted cell is a closed polyhedron, i.e. every edge of its face set is
# used exactly twice. It is checked TWICE, on purpose:
#   1. in-pipeline, via NINJA_CORE_INTEGRITY, so a stage that breaks it
#      is named even when a later stage happens to remove the culprit;
#   2. on the written polyMesh, via Benchmarks/check_integrity.py.
# A cell size floor is asserted too, so "remove everything" can never
# pass this gate.
#
# Usage: test_gatehanging.sh <ninjaMesher-exe> <caseDir>
set -eu

EXE="$1"
CASE_DIR="$2"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

OUT_DIR="$CASE_DIR/constant/polyMesh"
rm -rf "$OUT_DIR"

LOG="$(mktemp)"
trap 'rm -f "$LOG"' EXIT
NINJA_CORE_INTEGRITY=1 mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" "$CASE_DIR" | tee "$LOG"

status=0

# 1. In-pipeline contract at every stage that reports it.
STAGES=$(grep -c '^coreIntegrity ' "$LOG" || true)
if [ "$STAGES" -lt 1 ]; then
    echo "FAIL: no coreIntegrity lines -- NINJA_CORE_INTEGRITY did not report"
    status=1
fi
if grep '^coreIntegrity ' "$LOG" | grep -v 'nonClosedCells=0 *$' > /dev/null; then
    echo "FAIL: a pipeline stage emitted non-edge-closed cells:"
    grep '^coreIntegrity ' "$LOG" | grep -v 'nonClosedCells=0 *$'
    status=1
fi

# 2. Same contract on the written mesh.
if ! python3 "$REPO_ROOT/Benchmarks/check_integrity.py" "$OUT_DIR"; then
    echo "FAIL: check_integrity.py reported a structural mesh pathology"
    status=1
fi

# 3. Floor on the kept mesh: the fix removes exactly the cells that
# cannot be emitted as closed polyhedra (1 of 493 as measured, leaving
# 492), so a wholesale collapse can never satisfy the gate by deleting
# the offending geometry.
NCELLS=$(sed -n 's/.*wrote [0-9]* points, \([0-9]*\) cells.*/\1/p' "$LOG" | tail -1)
echo "keptCells=$NCELLS"
if [ -z "$NCELLS" ] || [ "$NCELLS" -lt 480 ]; then
    echo "FAIL: only $NCELLS cells kept (>= 480 expected; 492 measured at the fix)"
    status=1
fi

exit $status
