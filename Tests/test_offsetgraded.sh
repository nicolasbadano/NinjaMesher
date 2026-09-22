#!/bin/sh
# Gate for offset_graded: the layer GRADING contract itself.
#
# Every earlier offset gate would stay green if the graded march
# collapsed back to a single full-thickness step -- none of them looks
# at layer geometry. This gate does:
#
#   * 4 graded layers, expansionRatio 1.3, finalLayerThickness 0.1775 on the
#     offset_sphere geometry (clearance chosen so nothing drops):
#     droppedFaces == 0 and frozenPoints == 0,
#   * perStepPrismCells has exactly 4 entries and they are all equal
#     (every step marched the whole front -- a collapsed march would
#     print one entry),
#   * measured grading: layer_grading.py walks wall-normal prism stacks
#     and reports mean successive-layer ratio (must be within +-15% of
#     1.3) and thin-at-wall monotonicity (must hold on every stack),
#   * determinism: a second run is byte-identical,
#   * np=1 vs np=4 byte-identical (layers run in the rank-0 serial tail,
#     so this is cheap but is the assertion that keeps it that way),
#   * the removed ABSOLUTE-thickness dict form is refused, and the
#     finalLayerThickness its error message quotes reproduces this mesh,
#   * absorbVolFrac is refused outside [0, 1] -- a negative value used to
#     be read as "not supplied" and silently fell back to the plain-cut
#     threshold, 300x smaller than the default,
#   * "Mesh OK." + check_integrity.py.
#
# Usage: test_offsetgraded.sh <ninjaMesher-exe> <caseDir>
set -eu

EXE="$1"
CASE_DIR="$2"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

python3 "$SCRIPT_DIR/scripts/make_stls.py" sphere_at 0.5 0.5 0.5 0.25 3 "$CASE_DIR/sphere.stl"

OUT_DIR="$CASE_DIR/constant/polyMesh"
NLAYERS=4
RATIO=1.3

status=0

rm -rf "$OUT_DIR"
mpirun --bind-to none -np 1 "$EXE" "$CASE_DIR" | tee "$TMP/run1.log"
cp -r "$OUT_DIR" "$TMP/mesh_run1"

echo "--- --cut-stats ---"
mpirun --bind-to none -np 1 "$EXE" --cut-stats "$CASE_DIR"

# --- stats: nothing dropped, nothing frozen ---------------------------
DROPPED=$(grep -E "^droppedFaces" "$TMP/run1.log" | sed -E 's/.*= ([0-9]+).*/\1/' || true)
FROZEN=$(grep -E "^frozenPoints" "$TMP/run1.log" | sed -E 's/.*= ([0-9]+).*/\1/' || true)
echo "measured droppedFaces = ${DROPPED:-<missing>}  frozenPoints = ${FROZEN:-<missing>}"
if [ "${DROPPED:-x}" != "0" ]; then
    echo "FAIL: droppedFaces expected 0, got ${DROPPED:-<missing>}"
    status=1
fi
if [ "${FROZEN:-x}" != "0" ]; then
    echo "FAIL: frozenPoints expected 0, got ${FROZEN:-<missing>}"
    status=1
fi

# --- per-step prism cells: nLayers entries, all equal -----------------
PSPC=$(grep -E "^perStepPrismCells" "$TMP/run1.log" | sed -E 's/^perStepPrismCells = //' || true)
echo "measured perStepPrismCells = ${PSPC:-<missing>}"
if [ -z "${PSPC:-}" ]; then
    echo "FAIL: 'perStepPrismCells' line not found (per-step stats missing)"
    status=1
else
    OK=$(echo "$PSPC" | awk -v n="$NLAYERS" \
        '{ if (NF != n) { print "0"; exit } for (i = 2; i <= NF; ++i) if ($i != $1) { print "0"; exit } print "1" }')
    if [ "$OK" != "1" ]; then
        echo "FAIL: perStepPrismCells must have $NLAYERS equal entries, got '$PSPC'"
        status=1
    fi
fi

# --- grading geometry -------------------------------------------------
GRADE=$(python3 "$SCRIPT_DIR/scripts/layer_grading.py" "$OUT_DIR" sphere "$NLAYERS" 40) || {
    echo "FAIL: layer_grading.py could not walk the prism stacks"
    status=1
    GRADE=""
}
echo "--- layer grading ---"
echo "$GRADE"
if [ -n "$GRADE" ]; then
    MEAN_RATIO=$(echo "$GRADE" | awk '/^meanRatio/ { print $2 }')
    MONOTONE=$(echo "$GRADE" | awk '/^monotone/ { print $2 }')
    OK=$(awk -v a="$MEAN_RATIO" -v r="$RATIO" \
        'BEGIN { print (a >= 0.85 * r && a <= 1.15 * r) ? "1" : "0" }')
    if [ "$OK" != "1" ]; then
        echo "FAIL: mean layer ratio $MEAN_RATIO not within +-15% of $RATIO"
        status=1
    fi
    if [ "$MONOTONE" != "1" ]; then
        echo "FAIL: layer heights are not monotonically thin-at-wall on every stack"
        status=1
    fi
fi

# --- determinism (double run) ----------------------------------------
rm -rf "$OUT_DIR"
mpirun --bind-to none -np 1 "$EXE" "$CASE_DIR" >"$TMP/run2.log" 2>&1
cp -r "$OUT_DIR" "$TMP/mesh_run2"
if diff -r "$TMP/mesh_run1" "$TMP/mesh_run2" >/dev/null 2>&1; then
    echo "PASS: determinism (two np=1 runs byte-identical)"
else
    echo "FAIL: two np=1 runs produced DIFFERENT polyMesh output"
    status=1
fi

# --- np-invariance (np=1 vs np=4) ------------------------------------
rm -rf "$OUT_DIR"
mpirun --bind-to none -np 4 "$EXE" "$CASE_DIR" >"$TMP/run_np4.log" 2>&1
cp -r "$OUT_DIR" "$TMP/mesh_np4"
if diff -r "$TMP/mesh_run1" "$TMP/mesh_np4" >/dev/null 2>&1; then
    echo "PASS: np=1 vs np=4 byte-identical"
else
    echo "FAIL: np=1 vs np=4 polyMesh output DIFFERS"
    status=1
fi

# --- absolute-thickness dict is refused, and its advice is CORRECT -----
# The absolute `totalThickness` grammar was removed on 2026-09-18. Two
# things must hold, and the second is the one that matters: the dict is
# refused, AND the `finalLayerThickness` the refusal quotes reproduces
# the mesh the old dict asked for. An error message carrying a wrong
# number would be worse than no message, since the whole point of the
# change was to stop people hand-converting thicknesses.
ABS_CASE="$TMP/absolute"
mkdir -p "$ABS_CASE"
cp -r "$CASE_DIR/system" "$CASE_DIR/sphere.stl" "$ABS_CASE/"
cp "$CASE_DIR/system/ninjaMeshDict.absolute" "$ABS_CASE/system/ninjaMeshDict"
if mpirun --bind-to none -np 1 "$EXE" "$ABS_CASE" >"$TMP/abs.log" 2>&1; then
    echo "FAIL: absolute 'totalThickness' dict was ACCEPTED; it must be refused"
    status=1
elif ! grep -q "no longer accepted" "$TMP/abs.log"; then
    echo "FAIL: absolute dict was refused without the migration message"
    cat "$TMP/abs.log"
    status=1
else
    SUGG=$(sed -nE "s/.*finalLayerThickness ([0-9.eE+-]+);.*/\1/p" "$TMP/abs.log" | head -1)
    echo "measured suggested finalLayerThickness = ${SUGG:-<missing>}"
    if [ -z "${SUGG:-}" ]; then
        echo "FAIL: refusal did not quote a finalLayerThickness to use"
        status=1
    else
        sed -E "s/finalLayerThickness[[:space:]]+[0-9.eE+-]+;/finalLayerThickness $SUGG;/" \
            "$CASE_DIR/system/ninjaMeshDict" >"$ABS_CASE/system/ninjaMeshDict"
        if ! mpirun --bind-to none -np 1 "$EXE" "$ABS_CASE" >"$TMP/abs2.log" 2>&1; then
            echo "FAIL: the finalLayerThickness quoted by the refusal does not run"
            cat "$TMP/abs2.log"
            status=1
        else
            for KEY in wallArea residualMean residualMax frozenPoints droppedFaces perStepPrismCells; do
                A=$(grep -E "^$KEY " "$TMP/run1.log" || true)
                B=$(grep -E "^$KEY " "$TMP/abs2.log" || true)
                if [ "$A" != "$B" ]; then
                    echo "FAIL: migrated dict stat mismatch: '$A' vs '$B'"
                    status=1
                fi
            done
            echo "PASS: absolute dict refused, and the quoted fraction reproduces the mesh"
        fi
    fi
fi

# --- absorbVolFrac bounds ---------------------------------------------
# `absorbVolFrac` is a volume FRACTION, so it means nothing outside
# [0, 1]. The negative case is the one that matters: cutMesh reads
# `volFracThreshold >= 0.0 ? volFracThreshold : kSmallVolFrac`, so before
# this validation a negative value was treated as "not supplied" and
# silently fell back to the plain-cut 1e-3 -- 300x smaller than the 0.3
# default, changing which cells get absorbed, printing nothing. The dict
# must be REFUSED rather than quietly meshed with a different threshold.
AVF_CASE="$TMP/absorbvolfrac"
mkdir -p "$AVF_CASE"
cp -r "$CASE_DIR/system" "$CASE_DIR/sphere.stl" "$AVF_CASE/"

# Writes the case dict with `absorbVolFrac <value>;` inside `layers`.
# The key goes after the block's OPENING BRACE, which is on its own line
# here -- inserting after the `layers` line itself produces a dict that
# fails to parse for the wrong reason, and every case would "pass" the
# refusal checks while telling us nothing.
write_avf() {
    awk -v v="$1" '
        /^layers/ { print; inlayers = 1; next }
        inlayers && /^\{/ { print; print "    absorbVolFrac " v ";"; inlayers = 0; next }
        { print }
    ' "$CASE_DIR/system/ninjaMeshDict" >"$AVF_CASE/system/ninjaMeshDict"
}

for BAD in -0.5 1.5; do
    write_avf "$BAD"
    if mpirun --bind-to none -np 1 "$EXE" "$AVF_CASE" >"$TMP/avf.log" 2>&1; then
        echo "FAIL: absorbVolFrac $BAD was ACCEPTED; it is outside [0, 1]"
        status=1
    elif ! grep -q "absorbVolFrac' must be in \[0, 1\]" "$TMP/avf.log"; then
        echo "FAIL: absorbVolFrac $BAD refused, but not by the bounds check"
        cat "$TMP/avf.log"
        status=1
    else
        echo "PASS: absorbVolFrac $BAD refused with the bounds message"
    fi
done

# The in-range endpoints and a normal value must still run. 0 is legal
# and disables absorption (no volFrac is ever below it).
for GOOD in 0 0.3 1; do
    write_avf "$GOOD"
    if ! mpirun --bind-to none -np 1 "$EXE" "$AVF_CASE" >"$TMP/avf.log" 2>&1; then
        echo "FAIL: absorbVolFrac $GOOD was refused; it is inside [0, 1]"
        cat "$TMP/avf.log"
        status=1
    else
        echo "PASS: absorbVolFrac $GOOD accepted"
    fi
done

# --- checkMesh + integrity (on the np=4 output, still in place) -------
rm -rf "$OUT_DIR"
cp -r "$TMP/mesh_run1" "$OUT_DIR"

if ! command -v checkMesh >/dev/null 2>&1; then
    echo "FAIL: checkMesh not found on PATH (source OpenFOAM environment first)"
    exit 1
fi

LOG="$TMP/checkMesh.log"
if ! checkMesh -case "$CASE_DIR" >"$LOG" 2>&1; then
    echo "FAIL: checkMesh exited non-zero"
    cat "$LOG"
    exit 1
fi
if ! grep -q "Mesh OK\." "$LOG"; then
    echo "FAIL: 'Mesh OK.' verdict not found in checkMesh output"
    status=1
fi
if ! python3 "$SCRIPT_DIR/../Benchmarks/check_integrity.py" "$OUT_DIR"; then
    echo "FAIL: check_integrity.py reported a structural pathology"
    status=1
fi

if [ "$status" -ne 0 ]; then
    cat "$LOG"
fi

exit $status
