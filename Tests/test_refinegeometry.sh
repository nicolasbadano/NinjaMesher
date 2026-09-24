#!/bin/sh
# `refinementGeometry`: an STL that only sizes the grid.
#
# A surface refinement band may name an STL that is not a wall -- the
# pieces a snappy case refines at different levels, say, when the wall
# itself is their closed union. Three assertions:
#   1. a band on a refinement-only COPY of the wall STL gives a mesh
#      byte-identical to the same band on the wall itself -- the same
#      refinement, and no extra patch, cut or layer from the copy;
#   2. `layers` on a refinement-only STL is refused;
#   3. one file declared in both blocks is refused.
#
# Usage: test_refinegeometry.sh <ninjaMesher-exe> <caseDir>
set -eu

EXE="$1"
CASE_DIR="$2"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
status=0

mesh() {  # <dict-file> <outDir>
    mkdir -p "$2"
    cp -rL "$CASE_DIR/system" "$CASE_DIR/sphere.stl" "$CASE_DIR/sizer.stl" "$2/"
    cp "$2/system/$1" "$2/system/ninjaMeshDict"
    mpirun --bind-to none -np 1 "$EXE" "$2" >"$2/log" 2>&1
}

# --- 1. refinement-only copy == the wall itself ------------------------
if mesh ninjaMeshDict.wall "$TMP/wall" && mesh ninjaMeshDict.sizing "$TMP/sizing"; then
    if diff -rq "$TMP/wall/constant/polyMesh" "$TMP/sizing/constant/polyMesh" >/dev/null; then
        echo "PASS: a refinement-only surface refines exactly like the wall, and adds nothing else"
    else
        echo "FAIL: mesh with the refinement-only copy differs from the wall-refined mesh"
        diff -rq "$TMP/wall/constant/polyMesh" "$TMP/sizing/constant/polyMesh" || true
        status=1
    fi
else
    echo "FAIL: a dict that must mesh did not"
    status=1
fi

# --- 2./3. refusals ----------------------------------------------------
refused() {  # <dict-file> <expected message fragment> <label>
    if mesh "$1" "$TMP/$1"; then
        echo "FAIL: $3 was ACCEPTED"
        status=1
    elif ! grep -q "$2" "$TMP/$1/log"; then
        echo "FAIL: $3 refused, but not with the expected message"
        cat "$TMP/$1/log"
        status=1
    else
        echo "PASS: $3 refused"
    fi
}
refused ninjaMeshDict.layersOnSizing "is a 'refinementGeometry' surface" "layers on a refinement-only STL"
refused ninjaMeshDict.both "declared in both 'geometry' and 'refinementGeometry'" "one file in both blocks"

exit $status
