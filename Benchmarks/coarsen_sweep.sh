#!/bin/bash
# NinjaMesher coarsening / over-thick-layer sweep.
#
# WHY: the contract this sweep exists to police is "when the user asks
# for cells too big for the geometry, or layers too thick for the cells,
# the mesh must NEVER come out with quality problems -- it should simply
# stay away from the wall". Under-resolution is ACCEPTABLE (a feature
# vanishes, the mesh is valid and away from the wall); an invalid cell
# is NOT, at any cell size.
#
# WHAT IT DOES: for a base case it generates scratch copies with
#   * grid axis      -- base grid n divided by 1, 2, 4, 8 (refinement
#                       `distance`/`radius` scaled by the same factor,
#                       per the convention bm_layers_fishpassage_coarse's
#                       dict header documents), layers ON and OFF;
#   * thickness axis -- at the SHIPPED grid, layer finalLayerThickness /
#                       firstLayerThickness multiplied by 1, 2, 4, 8
#                       ("too big" also means layers thicker than the
#                       local cell).
# and runs ninjaMesher -> check_integrity.py -> checkMesh -allGeometry
# -> summary.py on each, then prints one table.
#
# The grid is coarsened by keeping domain.min and the exact cell size
# multiple: dx_new = factor * dx_base, n_new = ceil(n/factor),
# max_new = min + n_new*dx_new. So every coarsened grid plane is a base
# grid plane (factor is a power of two) and the domain still covers the
# original -- the sub-box/grid-aligned discipline this suite uses.
#
# Scratch cases are written as Benchmarks/cases/sweep_<case>_<tag> --
# SAME DIRECTORY DEPTH as the real cases on purpose, so the relative
# STL symlinks some cases use (bm_layers_fishpassage_coarse/fp.stl ->
# ../bm_layers_fishpassage/fp.stl) still resolve. `Benchmarks/cases/
# sweep_*` is gitignored. run_all.sh names its cases explicitly, so
# these are invisible to it.
#
# Usage:
#   coarsen_sweep.sh [options] <case> [<case> ...]
#     --factors "1 2 4 8"     grid divisors           (default "1 2 4 8")
#     --thickness "1 2 4 8"   layer-thickness factors (default "1 2 4 8")
#     --thickness-grid <f>    run the thickness axis at grid factor <f>
#                             instead of the shipped grid (default 1).
#                             The question that axis asks -- what happens
#                             when the requested layer stack is thicker
#                             than the cell it grows from -- is a RATIO,
#                             so it is asked just as well on a coarser
#                             grid, at a fraction of the wall clock, on
#                             a case whose shipped grid is too heavy to
#                             sweep four times.
#     --no-grid               skip the grid axis
#     --no-thickness          skip the thickness axis
#     --layers-off-only       grid axis without the layered variants
#     --keep                  keep scratch case dirs (default: keep;
#                             --clean removes them after the run)
#     --clean                 remove each scratch case after measuring
#     --exe <path>            ninjaMesher (default build/ninjaMesher)
#     --np <n>                MPI ranks (default 1)
#     --table-only            re-print the table from existing results
#
# Results: Benchmarks/results-sweep/<tag>.{toml,mesher.log,checkmesh.log,
# integrity.log}. Never touches Benchmarks/results/ or gates.toml.
set -eu

BM_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$BM_DIR/.." && pwd)"
EXE="$REPO_DIR/build/ninjaMesher"
RESULTS_DIR="$BM_DIR/results-sweep"
GATES="$BM_DIR/gates.toml"
FACTORS="1 2 4 8"
THICKS="1 2 4 8"
THICK_GRID=1
DO_GRID=1
DO_THICK=1
LAYERS_VARIANTS="off on"
CLEAN=0
NP=1
TABLE_ONLY=0

while [ $# -gt 0 ]; do
    case "$1" in
        --factors)        FACTORS="${2:?}"; shift 2 ;;
        --thickness)      THICKS="${2:?}"; shift 2 ;;
        --thickness-grid) THICK_GRID="${2:?}"; shift 2 ;;
        --no-grid)        DO_GRID=0; shift ;;
        --no-thickness)   DO_THICK=0; shift ;;
        --layers-off-only) LAYERS_VARIANTS="off"; shift ;;
        --keep)           CLEAN=0; shift ;;
        --clean)          CLEAN=1; shift ;;
        --exe)            EXE="${2:?}"; shift 2 ;;
        --np)             NP="${2:?}"; shift 2 ;;
        --table-only)     TABLE_ONLY=1; shift ;;
        --) shift; break ;;
        -*) echo "unknown option $1" >&2; exit 2 ;;
        *)  break ;;
    esac
done

[ $# -ge 1 ] || { sed -n '2,50p' "$0"; exit 2; }
CASES="$*"

mkdir -p "$RESULTS_DIR"

# ---------------------------------------------------------------- dict
# Emits a variant ninjaMeshDict on stdout.
#   make_dict <src_dict> <grid_factor> <thickness_factor> <layers on|off>
make_dict() {
    python3 - "$1" "$2" "$3" "$4" <<'PY'
import math, re, sys

src, gfac, tfac, layers = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
txt = open(src).read()


def block_span(text, name):
    """(start, end) of `name { ... }` -- brace-matched, comments-tolerant."""
    m = re.search(r"(?m)^\s*" + re.escape(name) + r"\s*(?:\n\s*)?\{", text)
    if not m:
        return None
    i = text.index("{", m.start())
    depth = 0
    for j in range(i, len(text)):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                return (m.start(), j + 1)
    raise SystemExit(f"unbalanced braces in block {name} of {src}")


def strip_comments(s):
    return re.sub(r"//[^\n]*", "", s)


# ---- domain: keep min, set dx_new = gfac*dx, n_new = ceil(n/gfac),
#      max_new = min + n_new*dx_new (>= original max, grid-plane aligned).
ds, de = block_span(txt, "domain")
dom = txt[ds:de]
dclean = strip_comments(dom)


def vec3(key, text):
    # line-anchored: an unanchored r"n\s*\(" would match the trailing
    # "n" of "min (" and read the wrong vector.
    m = re.search(r"(?m)^\s*" + re.escape(key) + r"\s*\(([^)]*)\)", text)
    if not m:
        raise SystemExit(f"no domain.{key} in {src}")
    return [float(x) for x in m.group(1).split()]


lo = vec3("min", dclean)
hi = vec3("max", dclean)
n = [int(x) for x in vec3("n", dclean)]
dx = [(hi[k] - lo[k]) / n[k] for k in range(3)]
nn = [max(1, math.ceil(n[k] / gfac)) for k in range(3)]
nhi = [lo[k] + nn[k] * dx[k] * gfac for k in range(3)]
newdom = dom
newdom = re.sub(r"(?m)^(\s*max\s*)\(([^)]*)\)",
                lambda m: "%s(%.10g %.10g %.10g)" % (m.group(1), *nhi), newdom, count=1)
newdom = re.sub(r"(?m)^(\s*n\s*)\(([^)]*)\)",
                lambda m: "%s(%d %d %d)" % (m.group(1), *nn), newdom, count=1)
txt = txt[:ds] + newdom + txt[de:]

# ---- refinement: scale the ABSOLUTE band lengths with the cell size.
#      `level` and box/sphere placement are untouched -- a box is a
#      region of space, not a length that scales with dx.
span = block_span(txt, "refinement")
if span and gfac != 1:
    s, e = span
    blk = re.sub(r"(?m)\b(distance|radius)(\s+)([0-9.eE+-]+)\s*;",
                 lambda m: "%s%s%.10g;" % (m.group(1), m.group(2), float(m.group(3)) * gfac),
                 txt[s:e])
    txt = txt[:s] + blk + txt[e:]

# ---- layers
span = block_span(txt, "layers")
if span:
    s, e = span
    if layers == "off":
        txt = txt[:s] + "// layers block removed by coarsen_sweep.sh (layers OFF variant)\n" + txt[e:]
    elif tfac != 1.0 or gfac != 1:
        # `finalLayerThickness` is a FRACTION of the local cell, so a
        # coarser grid would already give proportionally thicker layers
        # for free. Dividing by gfac holds the ABSOLUTE thickness
        # constant along the grid axis, which is what this study measures
        # and what its published rows were measured at. The thickness axis
        # then scales by tfac on top.
        fac = tfac / float(gfac)
        blk = re.sub(r"(?m)\b(finalLayerThickness)(\s+)([0-9.eE+-]+)\s*;",
                     lambda m: "%s%s%.10g;" % (m.group(1), m.group(2), float(m.group(3)) * fac),
                     txt[s:e])
        txt = txt[:s] + blk + txt[e:]

hdr = ("// GENERATED by Benchmarks/coarsen_sweep.sh from %s\n"
       "// grid factor %d (dx x%d), layer thickness factor %g, layers %s\n"
       "// DO NOT EDIT -- scratch case, regenerate with the script.\n") % (src, gfac, gfac, tfac, layers)
sys.stdout.write(hdr + txt)
PY
}

# ------------------------------------------------------------- one run
#   run_variant <base_case> <tag> <grid_factor> <thick_factor> <on|off>
run_variant() {
    base="$1"; tag="$2"; gfac="$3"; tfac="$4"; lay="$5"
    src_dir="$BM_DIR/cases/$base"
    dst_dir="$BM_DIR/cases/sweep_$tag"
    [ -d "$src_dir" ] || { echo "no such case: $src_dir" >&2; return 1; }

    rm -rf "$dst_dir"
    mkdir -p "$dst_dir/system"
    # -a keeps the relative STL symlinks as symlinks; the scratch dir is
    # at the same depth as the source, so they still resolve.
    find "$src_dir" -maxdepth 1 -mindepth 1 ! -name constant ! -name system ! -name 'case.foam' \
        -exec cp -a {} "$dst_dir/" \;
    for f in "$src_dir"/system/*; do
        [ -e "$f" ] || continue
        case "$(basename "$f")" in ninjaMeshDict) ;; *) cp -a "$f" "$dst_dir/system/" ;; esac
    done
    make_dict "$src_dir/system/ninjaMeshDict" "$gfac" "$tfac" "$lay" >"$dst_dir/system/ninjaMeshDict"

    mlog="$RESULTS_DIR/${tag}.mesher.log"
    clog="$RESULTS_DIR/${tag}.checkmesh.log"
    ilog="$RESULTS_DIR/${tag}.integrity.log"
    t0=$(date +%s.%N)
    rc=0
    mpirun --bind-to none -np "$NP" "$EXE" "$dst_dir" >"$mlog" 2>&1 || rc=$?
    t1=$(date +%s.%N)
    elapsed=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f", b-a}')

    integ="skip"
    if [ "$rc" -eq 0 ]; then
        if python3 "$BM_DIR/check_integrity.py" "$dst_dir/constant/polyMesh" >"$ilog" 2>&1; then
            integ="ok"
        else
            integ="FAIL"
        fi
        checkMesh -case "$dst_dir" -allGeometry >"$clog" 2>&1 || true
    else
        integ="mesher-rc=$rc"
        : >"$clog"
        echo "ninjaMesher exited $rc; checkMesh not run" >>"$clog"
    fi

    # summary.py writes results-sweep/<tag>.toml (the full checkMesh
    # quality record). No [<tag>] section exists in gates.toml, so no
    # gate binds here -- this sweep MEASURES, run_all.sh gates.
    python3 "$BM_DIR/summary.py" "$tag" "$clog" "$GATES" "$RESULTS_DIR" "$elapsed" "$rc" "$mlog" \
        >/dev/null 2>&1 || true
    {
        echo "sweep_base_case = \"$base\""
        echo "sweep_grid_factor = $gfac"
        echo "sweep_thickness_factor = $tfac"
        echo "sweep_layers = \"$lay\""
        echo "sweep_integrity = \"$integ\""
        echo "sweep_mesher_rc = $rc"
    } >>"$RESULTS_DIR/${tag}.toml"
    printf '  %-44s rc=%s integrity=%-5s t=%ss\n' "$tag" "$rc" "$integ" "$elapsed"
    [ "$CLEAN" -eq 0 ] || rm -rf "$dst_dir"
}

# ------------------------------------------------------------- the run
if [ "$TABLE_ONLY" -eq 0 ]; then
    if [ ! -x "$EXE" ]; then echo "FAIL: no ninjaMesher at $EXE" >&2; exit 1; fi
    command -v checkMesh >/dev/null 2>&1 || { echo "FAIL: checkMesh not on PATH" >&2; exit 1; }
    echo "=== NinjaMesher coarsening sweep ==="
    echo "results -> $RESULTS_DIR"
    for base in $CASES; do
        echo "--- $base"
        if [ "$DO_GRID" -eq 1 ]; then
            for g in $FACTORS; do
                for lay in $LAYERS_VARIANTS; do
                    run_variant "$base" "${base}_g${g}_lay${lay}" "$g" 1 "$lay" || true
                done
            done
        fi
        if [ "$DO_THICK" -eq 1 ]; then
            for t in $THICKS; do
                # t=1 at the thickness-axis grid is the same run as that
                # grid's layers-ON variant; skip the duplicate.
                if [ "$t" = "1" ] && [ "$DO_GRID" -eq 1 ]; then
                    case " $FACTORS " in *" $THICK_GRID "*) continue ;; esac
                fi
                run_variant "$base" "${base}_g${THICK_GRID}_t${t}" "$THICK_GRID" "$t" on || true
            done
        fi
    done
fi

# --------------------------------------------------------------- table
python3 - "$RESULTS_DIR" $CASES <<'PY'
import os, re, sys

rd, cases = sys.argv[1], sys.argv[2:]


def read_record(path):
    """Lenient `key = value` reader for summary.py's per-case record.

    NOT tomllib: summary.py emits its free-text `failed_checks` /
    `checkmesh_warnings` values unquoted, which no TOML parser accepts.
    That file is a report, not a config, and it is summary.py's to
    change -- this sweep only reads it, so it reads it leniently.
    """
    rec = {}
    for line in open(path, errors="replace"):
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        k, v = k.strip(), v.strip()
        if v[:1] == '"' and v[-1:] == '"':
            rec[k] = v[1:-1]
        elif v in ("true", "false"):
            rec[k] = v == "true"
        else:
            try:
                rec[k] = int(v)
            except ValueError:
                try:
                    rec[k] = float(v)
                except ValueError:
                    rec[k] = v
    return rec


def mlog_num(txt, key, cast=float):
    m = re.search(r"(?m)^\s*" + re.escape(key) + r"\s*=\s*([0-9.eE+-]+)\s*$", txt)
    return cast(m.group(1)) if m else None


def fmt(v, spec="%.4g"):
    return "-" if v is None else (spec % v if isinstance(v, float) else str(v))


ACCEPTED = ("Concave cells", "small interpolation weight", "small volume ratio",
            "concave angles between consecutive edges", "small determinant")

hdr = ["variant", "cells", "t[s]", "verdict", "nonWL errors", "skew", "nonOrth",
       "aspect", "volume", "wallArea", "infeasWallA", "sealedA", "offsetCutComp", "cutDropComp", "cutDropVol", "integ"]
rows = []
for tag in sorted(os.listdir(rd)):
    if not tag.endswith(".toml"):
        continue
    tag = tag[:-5]
    if not any(tag.startswith(c + "_g") for c in cases):
        continue
    rec = read_record(os.path.join(rd, tag + ".toml"))
    mtxt = ""
    p = os.path.join(rd, tag + ".mesher.log")
    if os.path.exists(p):
        mtxt = open(p, errors="replace").read()
    failed = [s.strip() for s in rec.get("failed_checks", "").split(";") if s.strip()]
    nonwl = [s for s in failed if not any(a in s for a in ACCEPTED)]
    if rec.get("sweep_mesher_rc", 0) != 0:
        verdict = "MESHER-RC%d" % rec["sweep_mesher_rc"]
    elif rec.get("mesh_ok"):
        verdict = "Mesh OK."
    elif rec.get("effective_mesh_ok"):
        verdict = "OK(accepted)"
    elif rec.get("checkmesh_ran"):
        verdict = "FAIL"
    else:
        verdict = "n/a"
    rows.append([
        tag.replace("bm_layers_", ""),
        fmt(rec.get("cells")),
        "%.1f" % rec.get("wall_clock_seconds", 0.0),
        verdict,
        (("%d: " % len(nonwl)) + nonwl[0][:38]) if nonwl else "0",
        fmt(rec.get("max_skewness")),
        fmt(rec.get("max_nonortho")),
        fmt(rec.get("max_aspect_ratio")),
        fmt(rec.get("total_volume"), "%.6g"),
        fmt(rec.get("wall_area"), "%.6g"),
        fmt(mlog_num(mtxt, "infeasibleWallArea")),
        fmt(mlog_num(mtxt, "smoothSealedRegionArea")),
        fmt(mlog_num(mtxt, "offsetCutDisconnectedComponents", int)),
        fmt(mlog_num(mtxt, "discardedComponents", int)),
        fmt(mlog_num(mtxt, "discardedVolume")),
        rec.get("sweep_integrity", "?"),
    ])

if not rows:
    print("(no results)")
    sys.exit(0)
w = [max(len(hdr[i]), max(len(r[i]) for r in rows)) for i in range(len(hdr))]
sep = "  "
print()
print(sep.join(h.ljust(w[i]) for i, h in enumerate(hdr)))
print(sep.join("-" * w[i] for i in range(len(hdr))))
last = None
for r in sorted(rows, key=lambda r: r[0]):
    base = r[0].split("_g")[0]
    if last is not None and base != last:
        print()
    last = base
    print(sep.join(r[i].ljust(w[i]) for i in range(len(hdr))))
print()
print("nonWL errors: checkMesh *** errors outside summary.py's accepted whitelist.")
print("A NON-ZERO value there, or integ != ok, is a CONTRACT BREACH: coarse input")
print("must give an under-resolved but VALID mesh, never an invalid cell.")
PY
