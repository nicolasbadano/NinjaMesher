#!/bin/sh
# Duplicate dict keys are REFUSED, and distinct keys still combine by max.
#
# WHY THIS EXISTS. `parseEntries` stored entries with `entries[key] = ...`,
# so a repeated key silently overwrote the earlier one while
# `order.push_back(key)` still recorded it twice -- the two halves of a
# DictEntry disagreed and the loss was invisible. MEASURED on this very
# case: two refinement regions both named `dup`, level 2 declared first
# and level 1 second, produced 7 672 cells, bit-identical to the level-1
# region ALONE, against the 60 648 the same two regions give under
# distinct names. The exposure is worse outside `refinement`:
# buildGeometryConfig iterates `order` for patch ordering, so a repeated
# geometry key would emit a patch twice from one STL entry.
#
# Two assertions, because refusing duplicates is only half the contract:
#   1. distinct names -> BOTH regions apply, combined by max, so the
#      result equals the finer region on its own (nothing is dropped);
#   2. a duplicate name -> hard parse error naming the key.
#
# Usage: test_dictdupkey.sh <ninjaMesher-exe> <caseDir>
set -eu

EXE="$1"
CASE_DIR="$2"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
status=0

WORK="$TMP/case"
mkdir -p "$WORK"
cp -rL "$CASE_DIR/system" "$CASE_DIR/sphere.stl" "$WORK/"

cells_of() {  # <dict-file> -> prints cell count, or nothing on failure
    cp "$WORK/system/$1" "$WORK/system/ninjaMeshDict"
    mpirun --bind-to none -np 1 "$EXE" "$WORK" 2>&1 |
        sed -nE 's/.*wrote [0-9]+ points, ([0-9]+) cells.*/\1/p' | head -1
}

# --- 1. distinct names: max-combination, nothing dropped ---------------
cp "$CASE_DIR/system/ninjaMeshDict" "$WORK/system/ninjaMeshDict.distinct"
BOTH=$(cells_of ninjaMeshDict.distinct)
ONLY3=$(cells_of ninjaMeshDict.single2)
echo "measured cells: distinct-names=${BOTH:-<none>} level2-alone=${ONLY3:-<none>}"
if [ -z "${BOTH:-}" ] || [ -z "${ONLY3:-}" ]; then
    echo "FAIL: a dict that must parse did not produce a mesh"
    status=1
elif [ "$BOTH" != "$ONLY3" ]; then
    echo "FAIL: distinct-name regions did not combine by max ($BOTH vs $ONLY3 cells)"
    status=1
else
    echo "PASS: distinct names both apply, combined by max"
fi

# --- 2. duplicate name: refused, and the message names the key ---------
cp "$WORK/system/ninjaMeshDict.dupkey" "$WORK/system/ninjaMeshDict"
if mpirun --bind-to none -np 1 "$EXE" "$WORK" >"$TMP/dup.log" 2>&1; then
    echo "FAIL: duplicate dict key was ACCEPTED; it must be a parse error"
    status=1
elif ! grep -q "duplicate key 'dup'" "$TMP/dup.log"; then
    echo "FAIL: duplicate key refused, but the message does not name the key"
    cat "$TMP/dup.log"
    status=1
else
    echo "PASS: duplicate dict key refused, message names the key"
fi

[ "$status" -eq 0 ] && echo "test_dictdupkey: ALL PASS"
exit "$status"
