#!/bin/sh
# NinjaMesher benchmark suite.
#
# Runs every case in Benchmarks/cases: ninjaMesher -> checkMesh
# -allGeometry -> summary.py (parse -> per-case results/<case>.toml +
# one summary-table row). Exits nonzero iff any case's gates.toml gate
# fails. The family-6 solver-depth cases additionally run a real
# solver on the mesh they just built and gate on its convergence
# (bm_solver_smoke: 3 icoFoam steps on the bm_sphere_coarse mesh;
# bm_solver_fp_interfoam: a decomposed interFoam VoF run on the layered
# fish passage).
#
# Slow-case convention: a case marked `slow = true` in
# gates.toml (family 5 scale cases, family 6 solver-depth cases) is
# SKIPPED by default -- this keeps the default suite under ~2 minutes
# wall-clock, which is what the `benchmarks` ctest gate runs. Pass
# --all to also run every slow case (scale cases additionally record
# wall-clock + peak RSS via `/usr/bin/time -v`).
#
# Usage: run_all.sh [--all] [--only <case>] [ninjaMesher-exe]
#   --only <case>  run just that one case (implies --all, so a
#                  slow=true case is not skipped). Everything else is
#                  silently passed over -- use it to iterate on a single
#                  heavy case without a full-suite wait. Only that case's
#                  files in results/ are cleared; a full run's records
#                  survive a one-case re-check.
#   Default exe: <repo>/build/ninjaMesher (built by the normal
#   `cmake --build build` flow; run that first if missing).
set -eu

RUN_ALL=0
ONLY_CASE=""
while [ $# -gt 0 ]; do
    case "$1" in
        --all)  RUN_ALL=1; shift ;;
        --only) ONLY_CASE="${2:?--only needs a case name}"; RUN_ALL=1; shift 2 ;;
        *)      break ;;
    esac
done

BM_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$BM_DIR/.." && pwd)"
EXE="${1:-$REPO_DIR/build/ninjaMesher}"
MAKE_STLS="$REPO_DIR/Tests/scripts/make_stls.py"
RESULTS_DIR="$BM_DIR/results"
GATES="$BM_DIR/gates.toml"

# Reads gates.toml for `slow = true` under [<case_name>]. Best-effort
# (no real TOML parser in POSIX sh): looks for the case's section and a
# `slow = true` line before the next `[section]`.
case_is_slow() {
    awk -v want="[$1]" '
        $0 == want { infile=1; next }
        /^\[/ { infile=0 }
        infile && $0 ~ /^slow[[:space:]]*=[[:space:]]*true/ { found=1 }
        END { exit(found ? 0 : 1) }
    ' "$GATES"
}

# With --only, every case but the named one is passed over.
case_selected() {
    [ -z "$ONLY_CASE" ] || [ "$ONLY_CASE" = "$1" ]
}

# Reads a numeric `<key> = <value>` from gates.toml's [<case>] section.
# Same best-effort sh/awk parsing as case_is_slow (summary.py does the
# real TOML read for the checkMesh gates; these solver-convergence
# thresholds are consumed here, in the runner, and are documented as
# such in gates.toml).
gate_num() {
    awk -v want="[$1]" -v key="$2" -v dflt="$3" '
        $0 == want { infile=1; next }
        /^\[/ { infile=0 }
        infile && $0 ~ ("^" key "[[:space:]]*=") {
            sub(/^[^=]*=[[:space:]]*/, "", $0); sub(/[[:space:]]*(#.*)?$/, "", $0)
            found=$0
        }
        END { print (found == "" ? dflt : found) }
    ' "$GATES"
}

if [ ! -x "$EXE" ]; then
    echo "FAIL: ninjaMesher executable not found/executable at $EXE (build first)"
    exit 1
fi
if ! command -v checkMesh >/dev/null 2>&1; then
    echo "FAIL: checkMesh not found on PATH (source the OpenFOAM environment first)"
    exit 1
fi

# A full run owns the whole directory and starts clean. A --only run does
# NOT: it clears just the files of the case it is about to rewrite, so
# iterating on one case cannot destroy the record of a long suite that ran
# before it (a full --all run is hours; losing it to a one-case re-check is
# the kind of accident that only has to happen once).
mkdir -p "$RESULTS_DIR"
if [ -z "$ONLY_CASE" ]; then
    rm -rf "$RESULTS_DIR"
    mkdir -p "$RESULTS_DIR"
else
    rm -f "$RESULTS_DIR/$ONLY_CASE".* 2>/dev/null || true
fi

OVERALL=0
SUITE_START=$(date +%s.%N)

run_mesher_case() {
    case_name="$1"
    time_it="${2:-0}"   # 1: wrap the mesher run in /usr/bin/time -v (scale cases)
    case_selected "$case_name" || return 0
    if case_is_slow "$case_name" && [ "$RUN_ALL" -eq 0 ]; then
        echo "$case_name              SKIP   slow=true (run with --all to include)"
        return 0
    fi
    case_dir="$BM_DIR/cases/$case_name"
    rm -rf "$case_dir/constant/polyMesh"
    t0=$(date +%s.%N)
    mesher_rc=0
    if [ "$time_it" -eq 1 ] && command -v /usr/bin/time >/dev/null 2>&1; then
        time_log="$RESULTS_DIR/${case_name}.time.log"
        /usr/bin/time -v mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" "$case_dir" \
            >"$RESULTS_DIR/${case_name}.mesher.log" 2>"$time_log" || mesher_rc=$?
        peak_rss_kb=$(grep "Maximum resident set size" "$time_log" | awk '{print $NF}')
        echo "  (peak RSS: ${peak_rss_kb:-?} KB, see $time_log)"
    else
        mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" "$case_dir" >"$RESULTS_DIR/${case_name}.mesher.log" 2>&1 || mesher_rc=$?
    fi
    log="$RESULTS_DIR/${case_name}.checkmesh.log"
    if [ "$mesher_rc" -eq 0 ]; then
        # Structural integrity gate: pathologies checkMesh
        # and volume/area gates provably miss (duplicate-position
        # points, index-level cell closure, degenerate/duplicate
        # faces). ALWAYS a hard failure — expected_fail does NOT waive
        # it (it waives quality verdicts, never structural validity).
        if ! python3 "$BM_DIR/check_integrity.py" "$case_dir/constant/polyMesh" \
                >"$RESULTS_DIR/${case_name}.integrity.log" 2>&1; then
            echo "$case_name              FAIL   mesh integrity (see results/${case_name}.integrity.log)"
            cat "$RESULTS_DIR/${case_name}.integrity.log"
            OVERALL=1
            return 1
        fi
        checkMesh -case "$case_dir" -allGeometry >"$log" 2>&1 || true
    else
        : >"$log"
        echo "ninjaMesher exited $mesher_rc; checkMesh not run" >>"$log"
    fi
    t1=$(date +%s.%N)
    elapsed=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f", b-a}')
    rc=0
    python3 "$BM_DIR/summary.py" "$case_name" "$log" "$GATES" "$RESULTS_DIR" "$elapsed" "$mesher_rc" \
        "$RESULTS_DIR/${case_name}.mesher.log" || rc=$?
    return $rc
}

# Generic solver-depth case runner.
#
#   run_solver_case <case> <solver> [np] [already_meshed]
#
#     np             1 (default) runs the solver serially; N > 1 runs
#                    `decomposePar -force` + `mpirun -np N <solver>
#                    -parallel`. Needed for bm_solver_fp_interfoam,
#                    whose 1.39 M-cell interFoam step costs ~25 s
#                    serially.
#     already_meshed 1 skips the mesher call, e.g. when the caller
#                    already ran run_mesher_case on this case to get the
#                    checkMesh gate. Default 0 (mesh here).
#
# The solver record goes to results/<case>.solver.toml, NOT
# results/<case>.toml -- the latter belongs to summary.py's checkMesh
# record, and a case can legitimately have both.
#
# PRE-STEP: if system/setFieldsDict exists, 0/alpha.water is restored
# from 0/alpha.water.orig and `setFields` is run, so the case is
# re-runnable in place.
#
# CONVERGENCE GATE (this is the point of family 6: a mesh that passes
# checkMesh must also RUN). FAIL if any of:
#   - the solver exits non-zero;
#   - "FOAM FATAL" appears in the log;
#   - the last reported time is short of the dict's `endTime`, or the log
#     carries no "End" marker -- i.e. the run did not finish its time loop.
#     This, not the step COUNT, is the "did the solver survive this mesh"
#     question: a better mesh lets the controller take bigger steps and so
#     reaches endTime in FEWER of them;
#   - fewer than `min_steps` time steps were reported. A low floor only,
#     against a solver that exits 0 having done essentially nothing; it is
#     NOT a convergence measure and must not be re-tightened around a
#     measured step count (see gates.toml);
#   - the max Courant number ever exceeds `max_courant`;
#   - more than `max_bounding_lines` "bounding ..." lines appear
#     (occasional bounding of omega/k is normal; a flood is divergence);
#   - |cumulative continuity error| exceeds `max_cum_cont_err`.
# Thresholds come from gates.toml's [<case>] section via gate_num; the
# defaults below are permissive so the existing icoFoam/simpleFoam cases
# keep their old exit-code + FOAM-FATAL-only behaviour.
# The `endTime` of a case's controlDict, or empty if it cannot be read.
# Anchored on the KEY, so `stopAt endTime;` (endTime as a value, no number
# after it) does not match.
case_end_time() {
    sed -nE 's/^[[:space:]]*endTime[[:space:]]+([0-9.eE+-]+)[[:space:]]*;.*/\1/p' \
        "$1/system/controlDict" 2>/dev/null | tail -1
}

run_solver_case() {
    case_name="$1"
    solver="$2"
    np="${3:-1}"
    already_meshed="${4:-0}"
    case_dir="$BM_DIR/cases/$case_name"
    case_selected "$case_name" || return 0
    if case_is_slow "$case_name" && [ "$RUN_ALL" -eq 0 ]; then
        echo "$case_name              SKIP   slow=true (run with --all to include)"
        return 0
    fi
    if ! command -v "$solver" >/dev/null 2>&1; then
        echo "$case_name              SKIP   $solver not on PATH"
        return 0
    fi
    if [ "$already_meshed" -eq 0 ]; then
        rm -rf "$case_dir/constant/polyMesh"
        mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" "$case_dir" >"$RESULTS_DIR/${case_name}.mesher.log" 2>&1
    fi
    rm -rf "$case_dir"/[1-9]* "$case_dir"/0.0* "$case_dir"/processor* 2>/dev/null || true

    log="$RESULTS_DIR/${case_name}.${solver}.log"
    prelog="$RESULTS_DIR/${case_name}.presteps.log"
    : >"$prelog"
    prestep_rc=0
    if [ -f "$case_dir/system/setFieldsDict" ]; then
        [ -f "$case_dir/0/alpha.water.orig" ] && cp "$case_dir/0/alpha.water.orig" "$case_dir/0/alpha.water"
        setFields -case "$case_dir" >>"$prelog" 2>&1 || prestep_rc=$?
    fi
    if [ "$np" -gt 1 ] && [ "$prestep_rc" -eq 0 ]; then
        # decomposePar has no -numberOfSubdomains option, so the count is
        # injected through a generated dict rather than by editing the
        # case's own (tracked) system/decomposeParDict.
        dpd="$RESULTS_DIR/${case_name}.decomposeParDict"
        {
            echo "FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }"
            echo "numberOfSubdomains $np;"
            echo "method scotch;"
        } >"$dpd"
        decomposePar -case "$case_dir" -force -decomposeParDict "$dpd" >>"$prelog" 2>&1 || prestep_rc=$?
    fi

    t0=$(date +%s.%N)
    solver_rc=0
    if [ "$prestep_rc" -ne 0 ]; then
        : >"$log"
        echo "pre-step failed (rc=$prestep_rc); $solver not run" >>"$log"
        solver_rc=$prestep_rc
    elif [ "$np" -gt 1 ]; then
        mpirun --bind-to none -np "$np" "$solver" -case "$case_dir" -parallel >"$log" 2>&1 || solver_rc=$?
    else
        "$solver" -case "$case_dir" >"$log" 2>&1 || solver_rc=$?
    fi
    t1=$(date +%s.%N)
    elapsed=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f", b-a}')

    min_steps=$(gate_num "$case_name" min_steps 0)
    max_courant=$(gate_num "$case_name" max_courant 1e30)
    max_bounding=$(gate_num "$case_name" max_bounding_lines 1000000)
    max_cont=$(gate_num "$case_name" max_cum_cont_err 1e30)

    fatal=0
    grep -q "FOAM FATAL" "$log" && fatal=1
    nsteps=$(grep -c "^Time = " "$log" || true)
    nbounding=$(grep -c "^bounding " "$log" || true)
    # Completion, not step count, is the real "the solver did not die"
    # signal. A mesh that lets the time-step controller take BIGGER steps
    # reaches endTime in FEWER of them -- an improvement that a `min_steps`
    # floor reads as a failure. So gate on arrival: the last reported time
    # must reach the dict's endTime, and the solver must have left its time
    # loop normally (OpenFOAM prints "End" only then). `min_steps` survives
    # as a low floor against a solver that exits successfully having done
    # essentially nothing.
    end_time=$(case_end_time "$case_dir")
    last_time=$(grep -E "^Time = " "$log" | tail -1 | sed -E 's/^Time = //' | tr -d '[:space:]')
    # Under adjustTimeStep the final time legitimately lands SHORT of
    # endTime -- OpenFOAM's own loop condition is
    #     Time::run() -> value() < (endTime_ - 0.5*deltaT_)
    # (OpenFOAM-v2406, src/OpenFOAM/db/Time/Time.C:785), so the loop exits
    # as soon as one more step would overshoot by less than half a step.
    # Measured on this suite: gate stops 1.13e-4 short with deltaT 2.2e-3,
    # descargador 2.1e-4 short with deltaT 5.5e-4 -- both inside half a
    # step, neither inside any ULP-scale slack. So transcribe the upstream
    # criterion instead of inventing a tolerance. `deltaT = ` is printed by
    # the adjustable-step solvers; the fixed-step ones keep the dict value.
    last_delta=$(grep -oE "^deltaT = [0-9.eE+-]+" "$log" | tail -1 | sed 's/^deltaT = //')
    [ -n "$last_delta" ] || last_delta=$(sed -nE \
        's/^[[:space:]]*deltaT[[:space:]]+([0-9.eE+-]+)[[:space:]]*;.*/\1/p' \
        "$case_dir/system/controlDict" 2>/dev/null | tail -1)
    saw_end_marker=0
    grep -qE "^End$" "$log" && saw_end_marker=1
    reached_end=0
    if [ -n "$end_time" ] && [ -n "$last_time" ]; then
        awk -v v="$last_time" -v e="$end_time" -v d="${last_delta:-0}" \
            'BEGIN{ exit !(v+0 >= e+0 - 0.5*(d+0)) }' && reached_end=1
    fi
    peak_courant=$(grep -o "Courant Number mean: [0-9.eE+-]* max: [0-9.eE+-]*" "$log" \
        | awk '{ if ($NF+0 > m) m = $NF+0 } END { printf "%g", m+0 }')
    peak_cont=$(grep -o "cumulative = [0-9.eE+-]*" "$log" \
        | awk '{ v = $NF+0; if (v < 0) v = -v; if (v > m) m = v } END { printf "%g", m+0 }')

    status="PASS"
    reason="ok"
    if [ "$solver_rc" -ne 0 ]; then
        status="FAIL"; reason="$solver exited $solver_rc"
    elif [ "$fatal" -ne 0 ]; then
        status="FAIL"; reason="FOAM FATAL in log"
    elif [ -n "$end_time" ] && [ "$reached_end" -ne 1 ]; then
        status="FAIL"; reason="stopped at Time = ${last_time:-<none>}, more than half a step short of endTime $end_time"
    elif [ "$saw_end_marker" -ne 1 ]; then
        status="FAIL"; reason="solver log has no \"End\" marker (time loop did not finish)"
    elif [ "$nsteps" -lt "$min_steps" ]; then
        status="FAIL"; reason="only $nsteps time steps (min_steps=$min_steps)"
    elif awk -v v="$peak_courant" -v m="$max_courant" 'BEGIN{exit !(v+0 > m+0)}'; then
        status="FAIL"; reason="max Courant $peak_courant > max_courant=$max_courant"
    elif [ "$nbounding" -gt "$max_bounding" ]; then
        status="FAIL"; reason="$nbounding bounding lines (max_bounding_lines=$max_bounding)"
    elif awk -v v="$peak_cont" -v m="$max_cont" 'BEGIN{exit !(v+0 > m+0)}'; then
        status="FAIL"; reason="cumulative continuity error $peak_cont > max_cum_cont_err=$max_cont"
    fi
    [ "$status" = "PASS" ] || OVERALL=1

    # last few residual lines (final-iteration convergence, for reporting)
    last_residuals=$(grep -E "Solving for (p|p_rgh|Ux|Uy|Uz|alpha)" "$log" | tail -6 | tr '\n' '|')
    printf '%-24s %-6s steps=%-4s t=%.2fs  (%s)\n' "$case_name" "$status" "$nsteps" "$elapsed" "$reason"
    {
        echo "case = \"$case_name\""
        echo "solver = \"$solver\""
        echo "solver_np = $np"
        echo "solver_exit_code = $solver_rc"
        echo "foam_fatal_found = $([ "$fatal" -ne 0 ] && echo true || echo false)"
        echo "time_steps_or_iters_reported = $nsteps"
        echo "end_time_requested = \"${end_time:-}\""
        echo "last_time_reported = \"${last_time:-}\""
        echo "last_delta_t = \"${last_delta:-}\""
        echo "reached_end_time = $([ "$reached_end" -eq 1 ] && echo true || echo false)"
        echo "saw_end_marker = $([ "$saw_end_marker" -eq 1 ] && echo true || echo false)"
        echo "max_courant_observed = $peak_courant"
        echo "bounding_lines = $nbounding"
        echo "max_cum_cont_err_observed = $peak_cont"
        echo "wall_clock_seconds = $elapsed"
        echo "gate_result = \"$status\""
        echo "gate_reason = \"$reason\""
        echo "last_residual_lines = \"$last_residuals\""
    } >"$RESULTS_DIR/${case_name}.solver.toml"
    rm -rf "$case_dir"/processor* 2>/dev/null || true
    [ "$status" = "PASS" ]
}

# Verifies `--refine-stats`' maxAdjacentLevelDiff == 1 (2:1 grading) for
# a case with a `refinement` block. The refinement family's third gate
# component (validity + volume come from run_mesher_case/gates.toml;
# this is the 2:1 grading check).
check_refine_stats() {
    case_name="$1"
    case_selected "$case_name" || return 0
    case_dir="$BM_DIR/cases/$case_name"
    out=$(mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" --refine-stats "$case_dir" 2>&1) || {
        echo "$case_name              FAIL   --refine-stats invocation failed"
        OVERALL=1
        return 1
    }
    diff=$(echo "$out" | grep -o 'maxAdjacentLevelDiff = [0-9]*' | grep -o '[0-9]*$')
    if [ "$diff" != "1" ]; then
        echo "$case_name              FAIL   maxAdjacentLevelDiff=$diff (expected 1)"
        OVERALL=1
        return 1
    fi
    echo "  (refine-stats: maxAdjacentLevelDiff=1, ok)"
    return 0
}

# Verifies --refine-stats reports a nonzero cell count at every level in
# `$2` (space-separated), e.g. "0 1 2" -- proves a distance band actually
# formed multiple levels, not just that grading is legal (
# bm_dist_sphere).
check_levels_nonzero() {
    case_name="$1"
    levels="$2"
    case_selected "$case_name" || return 0
    case_dir="$BM_DIR/cases/$case_name"
    out=$(mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" --refine-stats "$case_dir" 2>&1) || {
        echo "$case_name              FAIL   --refine-stats invocation failed"
        OVERALL=1
        return 1
    }
    ok=1
    for lvl in $levels; do
        count=$(echo "$out" | grep -o "$lvl:[0-9]*" | grep -o '[0-9]*$' | head -1)
        if [ -z "$count" ] || [ "$count" -eq 0 ] 2>/dev/null; then
            echo "$case_name              FAIL   level $lvl count is zero/missing in --refine-stats"
            ok=0
        fi
    done
    if [ "$ok" -ne 1 ]; then
        OVERALL=1
        return 1
    fi
    echo "  (refine-stats: nonzero counts at levels $levels, ok)"
    return 0
}

echo "=== NinjaMesher benchmark suite ==="
echo "results -> $RESULTS_DIR"
echo

# --- bm_sphere_coarse ---
python3 "$MAKE_STLS" sphere_at 0.5 0.5 0.5 0.3 3 "$BM_DIR/cases/bm_sphere_coarse/sphere.stl"
run_mesher_case bm_sphere_coarse || OVERALL=1

# --- bm_sphere_refined ---
python3 "$MAKE_STLS" sphere_at 0.5 0.5 0.5 0.3 3 "$BM_DIR/cases/bm_sphere_refined/sphere.stl"
run_mesher_case bm_sphere_refined || OVERALL=1

# --- bm_two_spheres ---
python3 "$MAKE_STLS" sphere_at 0.3 0.5 0.5 0.12 3 "$BM_DIR/cases/bm_two_spheres/sphere1.stl"
python3 "$MAKE_STLS" sphere_at 0.7 0.5 0.5 0.12 3 "$BM_DIR/cases/bm_two_spheres/sphere2.stl"
run_mesher_case bm_two_spheres || OVERALL=1

# --- bm_rotcube_fine ---
python3 "$MAKE_STLS" rotcube "$BM_DIR/cases/bm_rotcube_fine/rotcube.stl"
run_mesher_case bm_rotcube_fine || OVERALL=1

# --- bm_thin_plate --- (thickness 0.4 * cell size 0.1 = 0.04, deliberately
# OFF a grid-vertex plane -- centred at x=0.53, not x=0.50 -- so both
# plate faces fall strictly inside the SAME cell-edge span (0.5-0.6)
# instead of landing on the cell boundary itself; this is what actually
# exercises the documented multi-cut limit, see final report.)
python3 "$MAKE_STLS" box_at 0.51 -0.05 -0.05 0.55 1.05 1.05 "$BM_DIR/cases/bm_thin_plate/plate.stl"
run_mesher_case bm_thin_plate || OVERALL=1

# --- bm_diag_jump --- (pure refinement, no STL)
run_mesher_case bm_diag_jump || OVERALL=1

# --- bm_sphere_offcenter ---
python3 "$MAKE_STLS" sphere_at 0.4718 0.5093 0.5211 0.25 3 "$BM_DIR/cases/bm_sphere_offcenter/sphere.stl"
run_mesher_case bm_sphere_offcenter || OVERALL=1

# --- bm_cavity_sphere ---
python3 "$MAKE_STLS" sphere_at 0.5 0.5 0.5 0.3 3 "$BM_DIR/cases/bm_cavity_sphere/sphere.stl"
run_mesher_case bm_cavity_sphere || OVERALL=1

# ============================================================
# Family 1 -- adverse placements
# ============================================================

# --- bm_adv_gridplane --- (sphere tangent to a grid plane + box face
# exactly on another grid plane)
python3 "$MAKE_STLS" sphere_at 0.3 0.5 0.5 0.2 3 "$BM_DIR/cases/bm_adv_gridplane/sphere.stl"
python3 "$MAKE_STLS" box_at 0.6 0.2 0.2 0.8 0.4 0.4 "$BM_DIR/cases/bm_adv_gridplane/box.stl"
run_mesher_case bm_adv_gridplane || OVERALL=1

# --- bm_adv_gridpoint --- (box vertices bit-exact on grid points)
python3 "$MAKE_STLS" box_at 0.3 0.3 0.3 0.6 0.6 0.6 "$BM_DIR/cases/bm_adv_gridpoint/box.stl"
run_mesher_case bm_adv_gridpoint || OVERALL=1

# --- bm_adv_tangent --- (sphere tangent to the domain boundary from inside)
python3 "$MAKE_STLS" sphere_at 0.3 0.5 0.5 0.3 3 "$BM_DIR/cases/bm_adv_tangent/sphere.stl"
run_mesher_case bm_adv_tangent || OVERALL=1

# --- bm_adv_edgegraze --- (box edge lies along a grid edge line)
python3 "$MAKE_STLS" box_at 0.5 0.5 0.2 0.8 0.8 0.5 "$BM_DIR/cases/bm_adv_edgegraze/box.stl"
run_mesher_case bm_adv_edgegraze || OVERALL=1

# ============================================================
# Family 2 -- STL x domain boundary
# ============================================================

# --- bm_bd_halfsphere --- (sphere centre on xmin face)
python3 "$MAKE_STLS" sphere_at 0.0 0.5 0.5 0.25 3 "$BM_DIR/cases/bm_bd_halfsphere/sphere.stl"
run_mesher_case bm_bd_halfsphere || OVERALL=1

# --- bm_bd_protrude --- (box protruding through ymax)
python3 "$MAKE_STLS" box_at 0.4 0.8 0.4 0.6 1.2 0.6 "$BM_DIR/cases/bm_bd_protrude/box.stl"
run_mesher_case bm_bd_protrude || OVERALL=1

# --- bm_bd_corner --- (sphere at a domain corner, 1/8 inside)
python3 "$MAKE_STLS" sphere_at 0.0 0.0 0.0 0.3 3 "$BM_DIR/cases/bm_bd_corner/sphere.stl"
run_mesher_case bm_bd_corner || OVERALL=1

# ============================================================
# Family 3 -- internal flow
# ============================================================

# --- bm_int_pipe --- (capped cylinder through the domain, fluid = bore)
python3 "$MAKE_STLS" capped_cylinder 0.5 0.5 -0.1 1.1 0.2 24 "$BM_DIR/cases/bm_int_pipe/pipe.stl"
run_mesher_case bm_int_pipe || OVERALL=1

# --- bm_int_elbow --- (single watertight compound L-tube, fluid = bore)
python3 "$MAKE_STLS" elbow_tube 0.2 0.5 0.2 0.5 0.5 0.2 0.5 0.5 0.8 0.08 16 "$BM_DIR/cases/bm_int_elbow/elbow.stl"
run_mesher_case bm_int_elbow || OVERALL=1

# ============================================================
# Family 4 -- refinement x geometry
# ============================================================

# --- bm_ref_partial --- (refinement covers only half the sphere surface)
python3 "$MAKE_STLS" sphere_at 0.5 0.5 0.5 0.3 3 "$BM_DIR/cases/bm_ref_partial/sphere.stl"
run_mesher_case bm_ref_partial || OVERALL=1
check_refine_stats bm_ref_partial || true

# --- bm_ref_level3 --- (level-3 region, base 8^3 grid)
python3 "$MAKE_STLS" sphere_at 0.5 0.5 0.5 0.3 3 "$BM_DIR/cases/bm_ref_level3/sphere.stl"
run_mesher_case bm_ref_level3 || OVERALL=1
check_refine_stats bm_ref_level3 || true

# --- bm_ref_multi --- (three disjoint regions, three different levels)
python3 "$MAKE_STLS" sphere_at 0.5 0.5 0.5 0.3 3 "$BM_DIR/cases/bm_ref_multi/sphere.stl"
run_mesher_case bm_ref_multi || OVERALL=1
check_refine_stats bm_ref_multi || true

# ============================================================
# Distance-mode surface refinement (per-octant octree)
# ============================================================

# --- bm_dist_sphere --- (10^3, icosphere r=0.3, distance band -> levels 0/1/2)
python3 "$MAKE_STLS" sphere_at 0.5 0.5 0.5 0.3 3 "$BM_DIR/cases/bm_dist_sphere/sphere.stl"
run_mesher_case bm_dist_sphere || OVERALL=1
check_refine_stats bm_dist_sphere || true
check_levels_nonzero bm_dist_sphere "0 1 2" || true

# --- bm_dist_thinplate --- (THE MARQUEE CASE: distance refinement makes
# the invisible-at-h=0.1 thin plate meshable, see final report)
python3 "$MAKE_STLS" box_at 0.51 -0.05 -0.05 0.55 1.05 1.05 "$BM_DIR/cases/bm_dist_thinplate/plate.stl"
run_mesher_case bm_dist_thinplate || OVERALL=1
check_refine_stats bm_dist_thinplate || true

# ============================================================
# Fast within-distance predicate
# ============================================================

# --- bm_dist_bigsoup --- (20^3, icosphere r=0.3 subdivision 6 -> 81920
# triangles, distance band 0.08/level 2 -- CI proxy for the exact
# early-exit `anyTriangleWithin` predicate's large-soup path without a
# heavy real-CAD import; HARD wall-clock gate, see gates.toml. Fast set
# (not slow=true): measured well under 60s -- see final report.)
python3 "$MAKE_STLS" sphere_at 0.5 0.5 0.5 0.3 6 "$BM_DIR/cases/bm_dist_bigsoup/sphere.stl"
run_mesher_case bm_dist_bigsoup 1 || OVERALL=1
check_refine_stats bm_dist_bigsoup || true
check_levels_nonzero bm_dist_bigsoup "0 1 2" || true

# ============================================================
# Imports: real CAD cases
# ============================================================

# --- bm_layers_fishpassage --- (layered sibling,
# completes 10m37s / 8.87 GB np=1; slow=true, expected_fail with the
# cited checkMesh quality mechanisms -- integrity is the hard gate)
run_mesher_case bm_layers_fishpassage 1 || OVERALL=1

# --- bm_layers_fishpassage_coarse --- (guiding-principle case, first
# gated after the layer-quality gate: 0 pyramids, 0 bad face
# tets, max skew 4.73 gated below 5.0; ~3 min np=1; slow=true)
run_mesher_case bm_layers_fishpassage_coarse 1 || OVERALL=1

# --- bm_layers_gate --- (real production spillway CAD, 4 inflation layers per surface,
# translated from the production addLayersControls; see case dict.)
run_mesher_case bm_layers_gate 1 || OVERALL=1

# --- bm_solver_gate_interfoam --- (interFoam VoF on the layered gated
# spillway: the fourth family-6 depth case. Two gates on ONE mesh, like
# bm_solver_fp_interfoam: the whitelisted-verdict checkMesh gate
# (run_mesher_case, same mesh as bm_layers_gate, whose ninjaMeshDict this
# case copies) AND a real interFoam convergence gate (run_solver_case,
# thresholds in gates.toml). 1.20M cells; the solver runs decomposed and
# reaches endTime 0.2 in 73 steps / ~265 s at np=8. slow=true.)
run_mesher_case bm_solver_gate_interfoam 1 || OVERALL=1
run_solver_case bm_solver_gate_interfoam interFoam "${NINJA_SOLVER_NP:-8}" 1 || OVERALL=1

# --- bm_layers_wfp --- (WFP fish-passage mouth, real CAD, 6 layers
# r=1.8 on fixedWalls; slow=true)
run_mesher_case bm_layers_wfp 1 || OVERALL=1

# --- bm_solver_wfp_interfoam --- (interFoam VoF on the layered WFP
# fish-passage mouth: the third family-6 depth case. Two gates on ONE
# mesh, like bm_solver_fp_interfoam: the whitelisted-verdict checkMesh
# gate (run_mesher_case, same thresholds as bm_layers_wfp, whose
# ninjaMeshDict this case copies) AND a real interFoam convergence gate
# (run_solver_case, thresholds in gates.toml). 1.09M cells; the solver
# runs decomposed, 35 steps to endTime 0.04 in ~95 s at np=8. slow=true.)
run_mesher_case bm_solver_wfp_interfoam 1 || OVERALL=1
run_solver_case bm_solver_wfp_interfoam interFoam "${NINJA_SOLVER_NP:-8}" 1 || OVERALL=1

# --- bm_layers_descargador --- (layered sibling, 4 layers r=1.5 on
# fixedWalls -- snappy asks 8, see the case dict for the measurement
# that motivated the reduction; 1.69M cells; slow=true)
run_mesher_case bm_layers_descargador 1 || OVERALL=1

# --- bm_solver_descargador_interfoam --- (interFoam VoF on the layered
# dam bottom outlet: the second family-6 depth case, same shape as
# bm_solver_fp_interfoam. Two gates on ONE mesh: the whitelisted-verdict
# checkMesh gate (run_mesher_case, identical thresholds to
# bm_layers_descargador's, whose ninjaMeshDict this case copies) AND a
# real interFoam convergence gate (run_solver_case, thresholds in
# gates.toml). 1.67M cells, so the solver runs decomposed: ~5.5 s/step
# at np=8, ~55 steps to endTime 0.04 in ~5 min. slow=true.)
run_mesher_case bm_solver_descargador_interfoam 1 || OVERALL=1
run_solver_case bm_solver_descargador_interfoam interFoam "${NINJA_SOLVER_NP:-8}" 1 || OVERALL=1

# ============================================================
# Family 5 -- scale (slow=true; skipped unless --all)
# ============================================================

# --- bm_scale_100 --- (100^3 = 1M cells, sphere, no refinement)
python3 "$MAKE_STLS" sphere_at 0.5 0.5 0.5 0.3 3 "$BM_DIR/cases/bm_scale_100/sphere.stl"
run_mesher_case bm_scale_100 1 || OVERALL=1

# --- bm_scale_refined --- (40^3 base + level-2 shell)
python3 "$MAKE_STLS" sphere_at 0.5 0.5 0.5 0.3 3 "$BM_DIR/cases/bm_scale_refined/sphere.stl"
run_mesher_case bm_scale_refined 1 || OVERALL=1

# --- bm_dist_scale --- (40^3 base + distance band, slow=true;
# what the triangle binning must survive -- see final report for
# wall-clock + peak RSS)
python3 "$MAKE_STLS" sphere_at 0.5 0.5 0.5 0.3 3 "$BM_DIR/cases/bm_dist_scale/sphere.stl"
run_mesher_case bm_dist_scale 1 || OVERALL=1

# ============================================================
# Family 6 -- solver depth (slow=true; skipped unless --all)
# ============================================================

# --- bm_solver_simple --- (simpleFoam, 200 iters laminar, bm_sphere_refined's mesh)
python3 "$MAKE_STLS" sphere_at 0.5 0.5 0.5 0.3 3 "$BM_DIR/cases/bm_solver_simple/sphere.stl"
run_solver_case bm_solver_simple simpleFoam || OVERALL=1

# --- bm_solver_ico_refined --- (icoFoam, 50 steps, bm_ref_partial's refined+cut mesh)
python3 "$MAKE_STLS" sphere_at 0.5 0.5 0.5 0.3 3 "$BM_DIR/cases/bm_solver_ico_refined/sphere.stl"
run_solver_case bm_solver_ico_refined icoFoam || OVERALL=1

# --- bm_solver_smoke: reuses bm_sphere_coarse's mesh
# (same ninjaMeshDict) + real 0/ fields + fvSchemes/fvSolution, runs 3
# icoFoam time steps. Gate: exit 0 and no "FOAM FATAL" in the log.
# --- bm_solver_fp_interfoam --- (interFoam VoF on the layered fish
# passage: the family-6 depth case. Two gates on ONE mesh: the
# whitelisted-verdict checkMesh gate (run_mesher_case, like
# bm_layers_fishpassage_coarse) AND a real interFoam convergence gate
# (run_solver_case, thresholds in gates.toml). Solver runs decomposed --
# ~4 s/step at np=8 vs ~25 s serial. slow=true.)
run_mesher_case bm_solver_fp_interfoam 1 || OVERALL=1
run_solver_case bm_solver_fp_interfoam interFoam "${NINJA_SOLVER_NP:-8}" 1 || OVERALL=1

SMOKE_DIR="$BM_DIR/cases/bm_solver_smoke"
if ! case_selected bm_solver_smoke; then
    :
else
python3 "$MAKE_STLS" sphere_at 0.5 0.5 0.5 0.3 3 "$SMOKE_DIR/sphere.stl"
rm -rf "$SMOKE_DIR/constant/polyMesh"
mpirun --bind-to none -np "${NINJA_NP:-1}" "$EXE" "$SMOKE_DIR" >"$RESULTS_DIR/bm_solver_smoke.mesher.log" 2>&1
if ! command -v icoFoam >/dev/null 2>&1; then
    echo "bm_solver_smoke        SKIP   icoFoam not on PATH"
else
    rm -rf "$SMOKE_DIR"/[1-9]* "$SMOKE_DIR"/0.0* 2>/dev/null || true
    t0=$(date +%s.%N)
    icofoam_rc=0
    icoFoam -case "$SMOKE_DIR" >"$RESULTS_DIR/bm_solver_smoke.icofoam.log" 2>&1 || icofoam_rc=$?
    t1=$(date +%s.%N)
    elapsed=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f", b-a}')
    fatal=0
    grep -q "FOAM FATAL" "$RESULTS_DIR/bm_solver_smoke.icofoam.log" && fatal=1
    status="PASS"
    reason="ok"
    if [ "$icofoam_rc" -ne 0 ]; then
        status="FAIL"
        reason="icoFoam exited $icofoam_rc"
        OVERALL=1
    elif [ "$fatal" -ne 0 ]; then
        status="FAIL"
        reason="FOAM FATAL in log"
        OVERALL=1
    fi
    nsteps=$(grep -c "^Time = " "$RESULTS_DIR/bm_solver_smoke.icofoam.log" || true)
    printf 'bm_solver_smoke_icofoam %-6s steps=%-4s t=%.2fs  (%s)\n' "$status" "$nsteps" "$elapsed" "$reason"
    {
        echo "case = \"bm_solver_smoke_icofoam\""
        echo "icofoam_exit_code = $icofoam_rc"
        echo "foam_fatal_found = $([ "$fatal" -ne 0 ] && echo true || echo false)"
        echo "time_steps_reported = $nsteps"
        echo "wall_clock_seconds = $elapsed"
        echo "gate_result = \"$status\""
    } >"$RESULTS_DIR/bm_solver_smoke_icofoam.toml"
fi
fi

SUITE_END=$(date +%s.%N)
SUITE_ELAPSED=$(awk -v a="$SUITE_START" -v b="$SUITE_END" 'BEGIN{printf "%.2f", b-a}')

echo
echo "=== suite wall-clock: ${SUITE_ELAPSED}s ==="
if [ "$OVERALL" -eq 0 ]; then
    echo "=== ALL BENCHMARK GATES PASSED ==="
else
    echo "=== SOME BENCHMARK GATES FAILED (see reasons above / results/*.toml) ==="
fi
exit $OVERALL
