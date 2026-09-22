#!/bin/sh
# Memory-scaling gate: per-rank peak RSS at np=4
# must be WELL BELOW the np=1 peak on a scale case (bm_scale_100, 100^3
# = ~1M cells). Gated on the NON-ROOT ranks (< 50% of np=1; measured
# ~11%): rank 0 is exempt by decision -- it assembles the full mesh at
# the gather and writes the serial-assembled case (the accepted rank-0
# ceiling), so its peak is
# ~the np=1 peak by design. Also byte-compares the np=4 output against
# np=1 (the invariance bar applies here too).
#
# Usage: test_rss_scaling.sh <ninjaMesher-exe> <repoRoot>
set -eu

EXE="$1"
ROOT="$2"
CASE_DIR="$ROOT/Benchmarks/cases/bm_scale_100"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

python3 "$ROOT/Tests/scripts/make_stls.py" sphere_at 0.5 0.5 0.5 0.3 3 "$CASE_DIR/sphere.stl"

rm -rf "$CASE_DIR/constant/polyMesh"
NINJA_RSS=1 mpirun --bind-to none -np 1 "$EXE" "$CASE_DIR" >"$TMP/np1.log" 2>"$TMP/np1.err"
cp -r "$CASE_DIR/constant/polyMesh" "$TMP/np1"
RSS1=$(grep "maxrss" "$TMP/np1.err" | awk '{print $4}')

rm -rf "$CASE_DIR/constant/polyMesh"
NINJA_RSS=1 mpirun --bind-to none -np 4 "$EXE" "$CASE_DIR" >"$TMP/np4.log" 2>"$TMP/np4.err"

status=0
if ! diff -r "$TMP/np1" "$CASE_DIR/constant/polyMesh" >/dev/null 2>&1; then
    echo "FAIL: bm_scale_100 np=4 output differs from np=1"
    status=1
fi

echo "np=1 peak RSS: $RSS1 KB"
LIMIT=$((RSS1 / 2))
grep "maxrss" "$TMP/np4.err" | while read -r _ rank _ rss _; do
    echo "np=4 rank $rank peak RSS: $rss KB"
done
for line in $(grep "maxrss" "$TMP/np4.err" | awk '{print $2":"$4}'); do
    rank=${line%%:*}
    rss=${line##*:}
    if [ "$rank" = "0" ]; then
        continue # rank-0 gather/writer ceiling, exempt by decision (see header)
    fi
    if [ "$rss" -ge "$LIMIT" ]; then
        echo "FAIL: rank $rank peak RSS $rss KB not below 50% of np=1 ($LIMIT KB)"
        status=1
    fi
done

if [ "$status" = 0 ]; then
    echo "PASS: all non-root np=4 ranks below 50% of the np=1 peak RSS"
fi
exit $status
