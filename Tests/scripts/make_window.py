#!/usr/bin/env python3
"""Extract a small window case around a bounding box reported by
locate_errors.py.

Given a parent case dir + a box (min/max in real coordinates, assumed
already snapped to the parent's base-mesh lattice by locate_errors.py),
emits a new case dir Tests/cases/win_<name> with:
  - domain = the box, same dx as the parent (n recomputed from the box
    extent / dx, so cell size is unchanged)
  - the SAME real STL(s), as relative symlinks
  - the parent's `layers`/`refinement` blocks copied verbatim (they are
    absolute-coordinate rules; shrinking the domain to the window
    naturally clips their effect to what falls inside it -- no
    resampling needed)
  - all boundary patches forced to `type patch` (the window's outer
    faces are an artificial cut through the parent domain, never a
    real domain wall)
  - a locationInMesh chosen by the rule below

locationInMesh rule (CORRECTED: the original
"first corner that validates as fluid" rule only checked that the
candidate was NOT solid -- it never checked that the candidate was in
the SAME fluid connected component the parent keeps. A window box can
straddle several disconnected fluid regions (e.g. two separate
channels either side of a wall, or fluid pockets split by a thin
STL feature); the first-corner-that-validates rule could -- and, on
`win_fish_coarse_pyramids`/`win_fish_coarse_site1`, DID -- silently
pick the WRONG region, keeping mesh on the opposite side of the STL
from what the parent case keeps. Fixed rule, in priority order:

  1. `--location x y z` (explicit CLI override): used verbatim, no
     probing at all. For manual control when the rule below still
     gets it wrong, or for reproducing a known-good point exactly.
  2. If the PARENT's own `locationInMesh` falls strictly inside the
     window box (pulled 1% of the box diagonal in from each face, same
     margin used for every other probe candidate below): use it
     directly. This is the highest-confidence choice because it is
     the exact point the parent case itself already validated as
     "the side we keep" -- no ambiguity possible.
  3. Otherwise, build a candidate set: the 8 corners of the box, the 6
     face centers, and (if `--defect-centroid` is given) the defect
     centroid itself clamped into the box plus its reflection through
     the box center clamped into the box -- all pulled 1% of the box
     diagonal inward. For EVERY candidate, run the real classifier
     (`ninjaMesher --cut-stats`) on a throwaway copy of the dict and
     keep the ones that succeed (rc=0). Among the successful
     candidates, the WINDOW IS RE-CUT ONCE PER CANDIDATE (seconds-
     scale on these window-sized boxes) and the candidate whose run
     reports the LARGEST `totalVolume` (equivalently `keptCells`, same
     dx so proportional) is selected -- this is the fluid component
     that occupies the most of the box, which for a window carved out
     of a real, mostly-single-component domain is overwhelmingly the
     side the parent case actually keeps. Ties (within 1e-9 relative)
     are broken by whichever candidate's direction from the box center
     is closer (larger dot product of unit vectors) to the parent's
     own `locationInMesh` direction from the box center, when the
     parent dict has one; failing that, by the original
     defect-centroid-distance ordering; failing that, by candidate
     index (deterministic either way).

Every evaluated candidate's `keptCells`/`totalVolume`/
`disconnectedCellsDropped`/`discardedComponents` (from `--cut-stats`)
is printed to stdout, so the selection is auditable in numbers, not
just by eyeballing the mesh afterward -- `discardedComponents` > 0 on
the FINAL choice is itself a signal the box spans multiple components
even after picking the largest one.

Usage:
  make_window.py <parentCaseDir> <name> \
      --min x y z --max x y z \
      [--defect-centroid x y z] [--location x y z] \
      [--exe path/to/ninjaMesher]

Emits Tests/cases/win_<name>.
"""
import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def read_text(p):
    return Path(p).read_text()


def parse_location(dict_text):
    m = re.search(r"locationInMesh\s*\(([^)]*)\)\s*;", dict_text)
    if not m:
        return None
    return tuple(float(x) for x in m.group(1).split())


def parse_domain(dict_text):
    m = re.search(r"domain\s*\{([^}]*)\}", dict_text, re.S)
    block = m.group(1)
    def vec(name):
        mm = re.search(name + r"\s*\(([^)]*)\)", block)
        return tuple(float(x) for x in mm.group(1).split())
    dmin = vec("min")
    dmax = vec("max")
    n = tuple(int(x) for x in re.search(r"(?<![a-zA-Z])n\s*\(([^)]*)\)", block).group(1).split())
    return dmin, dmax, n


def build_dict(parent_text, box_min, box_max, n, patches_block, extra_header):
    domain_block = (
        "domain\n{\n"
        f"    min ({box_min[0]:.10g} {box_min[1]:.10g} {box_min[2]:.10g});\n"
        f"    max ({box_max[0]:.10g} {box_max[1]:.10g} {box_max[2]:.10g});\n"
        f"    n   ({n[0]} {n[1]} {n[2]});\n"
        "}\n"
    )
    text = re.sub(r"domain\s*\{[^}]*\}", domain_block, parent_text, count=1, flags=re.S)
    # patches' block has NESTED braces (each patch entry is itself
    # "side { type ...; name ...; }"), so a [^}]* class (which cannot
    # span any '}' at all) stops at the first entry's closing brace.
    # Match up to the first STANDALONE "}" line instead (the project's
    # dict convention: every top-level block closes with '}' alone on
    # its own line).
    text = re.sub(r"patches\s*\{.*?\n\}\n", patches_block, text, count=1, flags=re.S)
    text = re.sub(r"locationInMesh\s*\([^)]*\)\s*;", "__LOCATION_PLACEHOLDER__", text, count=1)
    return extra_header + text


def patches_all_patch(parent_text):
    m = re.search(r"patches\s*\{(.*?)\n\}\n", parent_text, re.S)
    block = m.group(1)
    entries = re.findall(r"(\w+)\s*\{\s*type\s+\w+\s*;\s*name\s+(\w+)\s*;\s*\}", block)
    lines = ["patches\n{\n"]
    for side, name in entries:
        lines.append(f"    {side} {{ type patch; name {name}; }}\n")
    lines.append("}\n")
    return "".join(lines)


def try_location(exe, case_dir, template_text, candidate):
    """Write `candidate` into the dict and run --cut-stats. Returns
    (ok, text2, stats) where stats is a dict of the parsed cut-stats
    numbers on success (None on failure/timeout)."""
    dict_path = case_dir / "system" / "ninjaMeshDict"
    text2 = template_text.replace(
        "__LOCATION_PLACEHOLDER__",
        f"locationInMesh ({candidate[0]:.10g} {candidate[1]:.10g} {candidate[2]:.10g});")
    dict_path.write_text(text2)
    try:
        r = subprocess.run([str(exe), "--cut-stats", str(case_dir)],
                            capture_output=True, text=True, timeout=120)
        if r.returncode != 0:
            return False, text2, None
        stats = {}
        for key in ("keptCells", "totalVolume", "disconnectedCellsDropped",
                    "discardedComponents"):
            m = re.search(key + r"\s*=\s*([-0-9.eE+]+)", r.stdout)
            stats[key] = float(m.group(1)) if m else None
        return True, text2, stats
    except Exception:
        return False, text2, None


def point_in_box(p, box_min, box_max, margin_frac=0.01):
    """Strictly inside the box, with the same 1%-of-diagonal margin
    used for every probe candidate (so a point sitting exactly on a
    face boundary is not treated as 'in the box')."""
    if p is None:
        return False
    for i in range(3):
        d = box_max[i] - box_min[i]
        lo = box_min[i] + margin_frac * d
        hi = box_max[i] - margin_frac * d
        if not (lo <= p[i] <= hi):
            return False
    return True


def gen_candidates(box_min, box_max, defect_centroid, pull=0.01):
    """8 corners + 6 face centers + (if given) the defect centroid and
    its reflection through the box center, all clamped/pulled 1% of
    the box diagonal inside the box. Deterministic order (used as the
    final tie-break when nothing else distinguishes two candidates)."""
    diag = tuple(box_max[i] - box_min[i] for i in range(3))
    lo = tuple(box_min[i] + pull * diag[i] for i in range(3))
    hi = tuple(box_max[i] - pull * diag[i] for i in range(3))
    mid = tuple(0.5 * (box_min[i] + box_max[i]) for i in range(3))

    candidates = []
    for cx in (lo[0], hi[0]):
        for cy in (lo[1], hi[1]):
            for cz in (lo[2], hi[2]):
                candidates.append((cx, cy, cz))
    candidates += [
        (lo[0], mid[1], mid[2]), (hi[0], mid[1], mid[2]),
        (mid[0], lo[1], mid[2]), (mid[0], hi[1], mid[2]),
        (mid[0], mid[1], lo[2]), (mid[0], mid[1], hi[2]),
    ]
    if defect_centroid:
        dc = defect_centroid
        clamped = tuple(min(max(dc[i], lo[i]), hi[i]) for i in range(3))
        reflected = tuple(min(max(2 * mid[i] - dc[i], lo[i]), hi[i]) for i in range(3))
        candidates.append(clamped)
        candidates.append(reflected)

    # de-duplicate (clamped/reflected can coincide with a corner or a
    # face center) while preserving first-seen order.
    seen = set()
    unique = []
    for c in candidates:
        key = tuple(round(v, 9) for v in c)
        if key not in seen:
            seen.add(key)
            unique.append(c)
    return unique


def choose_location(exe, win_dir, dict_text, box_min, box_max, defect_centroid,
                     parent_location):
    """Implements the corrected locationInMesh rule -- see the module
    docstring. Returns (chosen_point, final_text, audit_rows) where
    audit_rows is a list of (candidate, ok, stats) for every point
    evaluated (empty if the parent's own location was reused
    directly)."""
    mid = tuple(0.5 * (box_min[i] + box_max[i]) for i in range(3))

    if point_in_box(parent_location, box_min, box_max):
        ok, text2, stats = try_location(exe, win_dir, dict_text, parent_location)
        if ok:
            print(f"  parent locationInMesh {parent_location} falls inside the "
                  f"window box and validates as fluid -- reused directly "
                  f"(highest-priority rule). stats={stats}")
            return parent_location, text2, [(parent_location, ok, stats)]
        print(f"  WARNING: parent locationInMesh {parent_location} is inside "
              f"the window box but does NOT validate (--cut-stats failed) -- "
              f"falling back to the candidate search.", file=sys.stderr)

    candidates = gen_candidates(box_min, box_max, defect_centroid)
    audit_rows = []
    results = []  # (candidate, text2, stats)
    for c in candidates:
        ok, text2, stats = try_location(exe, win_dir, dict_text, c)
        audit_rows.append((c, ok, stats))
        if ok:
            results.append((c, text2, stats))

    for c, ok, stats in audit_rows:
        print(f"  candidate {c}: ok={ok} stats={stats}")

    if not results:
        return None, dict_text, audit_rows

    def parent_dir_score(c):
        if parent_location is None:
            return 0.0
        v1 = tuple(c[i] - mid[i] for i in range(3))
        v2 = tuple(parent_location[i] - mid[i] for i in range(3))
        n1 = sum(x * x for x in v1) ** 0.5
        n2 = sum(x * x for x in v2) ** 0.5
        if n1 == 0 or n2 == 0:
            return 0.0
        return sum(v1[i] * v2[i] for i in range(3)) / (n1 * n2)

    def centroid_dist_score(c):
        if not defect_centroid:
            return 0.0
        return sum((c[i] - defect_centroid[i]) ** 2 for i in range(3))

    best_vol = max(r[2]["totalVolume"] for r in results)
    tol = max(1e-9 * abs(best_vol), 1e-12)
    tied = [r for r in results if abs(r[2]["totalVolume"] - best_vol) <= tol]
    tied.sort(key=lambda r: (-parent_dir_score(r[0]), -centroid_dist_score(r[0])))
    chosen_c, chosen_text, chosen_stats = tied[0]
    print(f"  SELECTED {chosen_c}: largest totalVolume={chosen_stats['totalVolume']} "
          f"among {len(results)} valid candidates "
          f"(discardedComponents={chosen_stats['discardedComponents']}, "
          f"disconnectedCellsDropped={chosen_stats['disconnectedCellsDropped']})")
    return chosen_c, chosen_text, audit_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("parentCaseDir")
    ap.add_argument("name")
    ap.add_argument("--min", nargs=3, type=float, required=True)
    ap.add_argument("--max", nargs=3, type=float, required=True)
    ap.add_argument("--defect-centroid", nargs=3, type=float, default=None)
    ap.add_argument("--location", nargs=3, type=float, default=None,
                    help="explicit locationInMesh override (skips all probing; "
                         "highest priority, see the corrected rule in the module "
                         "docstring)")
    ap.add_argument("--exe", default=str(REPO / "build" / "ninjaMesher"))
    ap.add_argument("--header", default=[], action="append",
                    help="extra header comment line to prepend (repeatable; "
                         "'// ' prefix and trailing newline added if missing)")
    args = ap.parse_args()

    parent_dir = Path(args.parentCaseDir).resolve()
    win_name = f"win_{args.name}"
    win_dir = REPO / "Tests" / "cases" / win_name

    parent_dict_text = read_text(parent_dir / "system" / "ninjaMeshDict")
    dmin, dmax, n = parse_domain(parent_dict_text)
    dx = tuple((dmax[i] - dmin[i]) / n[i] for i in range(3))

    box_min = tuple(args.min)
    box_max = tuple(args.max)
    win_n = tuple(max(1, round((box_max[i] - box_min[i]) / dx[i])) for i in range(3))

    if win_dir.exists():
        shutil.rmtree(win_dir)
    (win_dir / "system").mkdir(parents=True)
    (win_dir / "constant").mkdir(parents=True)

    # Copy non-geometry system files verbatim (controlDict etc, if present).
    for fname in ("controlDict", "fvSchemes", "fvSolution"):
        src = parent_dir / "system" / fname
        if src.exists():
            shutil.copy(src, win_dir / "system" / fname)
    case_foam = parent_dir / "case.foam"
    if case_foam.exists():
        (win_dir / "case.foam").touch()

    # Symlink every geometry STL referenced in the parent dict (the
    # convention: relative symlinks, not copies).
    stl_names = re.findall(r"(\S+\.stl)\s*\n?\s*\{", parent_dict_text)
    for stl in set(stl_names):
        src = parent_dir / stl
        if src.exists():
            rel = Path("..") / parent_dir.relative_to(REPO) / stl
            # relative path from win_dir to parent stl
            relpath = Path(
                __import__("os").path.relpath(src, win_dir)
            )
            (win_dir / stl).symlink_to(relpath)

    patches_block = patches_all_patch(parent_dict_text)
    extra = "".join(
        (line if line.startswith("//") else "// " + line).rstrip("\n") + "\n"
        for line in args.header
    )
    header = (
        "// Auto-generated by Tests/scripts/make_window.py.\n"
        f"// Parent case: {parent_dir.relative_to(REPO)}\n"
        f"// Extraction box (same dx as parent): min{box_min} max{box_max}\n"
        + extra
    )
    dict_text = build_dict(parent_dict_text, box_min, box_max, win_n, patches_block, header)
    (win_dir / "system" / "ninjaMeshDict").write_text(dict_text)

    # locationInMesh: corrected rule (see module docstring) -- explicit
    # --location override, else parent's own locationInMesh if it falls
    # inside the box, else a same-fluid-component candidate search
    # picking the largest kept volume.
    exe = Path(args.exe)
    template_text = dict_text
    parent_location = parse_location(parent_dict_text)
    diag = tuple(box_max[i] - box_min[i] for i in range(3))

    if args.location:
        chosen = tuple(args.location)
        ok, _, stats = try_location(exe, win_dir, template_text, chosen)
        print(f"  --location override {chosen}: ok={ok} stats={stats}")
        if not ok:
            print(f"WARNING: --location override {chosen} does NOT validate "
                  f"via --cut-stats (used anyway, per explicit manual "
                  f"override).", file=sys.stderr)
    else:
        chosen, _, audit_rows = choose_location(
            exe, win_dir, template_text, box_min, box_max,
            args.defect_centroid, parent_location)
        if chosen is None:
            print("WARNING: no candidate locationInMesh validated via "
                  "--cut-stats; using box center as a last resort (may be "
                  "inside solid).", file=sys.stderr)
            chosen = tuple(box_min[i] + 0.5 * diag[i] for i in range(3))
            try_location(exe, win_dir, template_text, chosen)

    # Whichever branch ran above may have left the on-disk dict pointing
    # at the LAST candidate it happened to try (the probing loop writes
    # every candidate in turn) rather than the SELECTED one -- pin it
    # down explicitly here so the emitted case always matches `chosen`.
    final_text = template_text.replace(
        "__LOCATION_PLACEHOLDER__",
        f"locationInMesh ({chosen[0]:.10g} {chosen[1]:.10g} {chosen[2]:.10g});")
    (win_dir / "system" / "ninjaMeshDict").write_text(final_text)

    print(f"{win_name}: locationInMesh = {chosen}, n = {win_n}")


if __name__ == "__main__":
    main()
