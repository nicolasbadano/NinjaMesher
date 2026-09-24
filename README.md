<p align="center">
  <img src="ninjaMesher.png" alt="NinjaMesher" width="320">
</p>

# NinjaMesher

NinjaMesher is a hex-dominant Cartesian cut-cell mesher for CFD, producing
[OpenFOAM](https://www.openfoam.com) `polyMesh` output directly. It takes a
background hex grid, refines it near geometry, cuts it against STL surfaces,
and extrudes prismatic inflation layers on wall patches — prioritizing mesh
quality (a clean `checkMesh` verdict) over exact geometric fidelity.

## Design: quality over exactness

NinjaMesher is built on a deliberate inversion of the usual meshing priority:
**cell quality is the constraint; geometric fidelity is the residual.** The
error from approximating the wall is bounded by the local cell size — visible,
and reducible with local refinement. The error from bad cell shapes (skew,
non-orthogonality, slivers) silently degrades every gradient and flux the
solver computes. So instead of snapping cells onto every ripple of the STL,
NinjaMesher cuts a perfect hex grid and keeps every cell solver-friendly by
construction: convex cut cells, one cut face per cell, bounded 2:1 refinement
grading.

This buys clean `checkMesh` verdicts on real CAD, and it costs exactness in
disclosed, bounded ways:

- **The wall is approximated, not reproduced.** Cut positions snap to the
  half-grid (grid vertices and edge midpoints), so the wall can deviate from
  the STL by up to a quarter of the local cell size. The remedy for tighter
  walls is local refinement, never a more exact cutter.
- **Sub-cell features do not exist.** One cut face per cell means a thin
  plate or gap smaller than a cell is not represented at all — not
  approximated coarsely, absent. This is what choosing a cell size means,
  and it is the same lever as everywhere else: resolve the features you
  care about with refinement, and they appear.
- **Layer landing is smoothed.** With `smoothRadius > 0` the inflation
  layers land on a smoothed version of the surface; sharp creases are
  rounded at the smoothing-radius scale, and features comparable to that
  radius can shrink or seal. What smoothing changes is measured and printed
  (`smoothSealedRegion*`, `smoothMoved*` statistics) — loud, never silent.
- **Layers degrade gracefully, not silently.** Where full-depth layers are
  geometrically infeasible, faces get locally fewer/thicker layers or none,
  and the affected wall area is reported (`infeasibleWallArea`,
  `droppedFaces`).

The default configuration (half-grid cut, planarized layer tops, smoothing
radius 1.0) was selected by sweeping the alternatives over the benchmark
suite and keeping the combination with the best checkMesh verdicts, layer
completeness, and wall-clock cost.

## Dependencies

- A POSIX platform. Linux is what this is developed and tested on; the
  code uses `<sys/resource.h>` for memory reporting, so it does not build
  on Windows without a POSIX layer.
- CMake ≥ 3.16, a C++17 compiler
- MPI (required at build time; single-rank runs are fine)
- OpenMP (optional, used if found)
- OpenFOAM — **not needed to build or run the mesher, mandatory to run the
  test suite.** Meshing a case requires none of it. But most `ctest` gates
  and the whole benchmark suite validate their output with `checkMesh`, and
  the solver gates run `interFoam` / `simpleFoam` / `icoFoam`, so they fail
  rather than skip if the OpenFOAM environment is not sourced. Developed
  and measured against OpenFOAM v2406 (openfoam.com).
- Python 3 (standard library only; test/benchmark harness scripts)

## Building

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

This produces `build/ninjaMesher` plus three unit-test binaries.

## Running

A *case* is a directory in OpenFOAM layout containing `system/ninjaMeshDict`
and the referenced STL files. The mesh is written to `<caseDir>/constant/polyMesh`.

The STL input is assumed watertight, non-self-intersecting and free of
coincident or overlapping bodies; it is not validated. Run `surfaceCheck` /
`surfaceClean` on anything that has not been through that already.

```sh
build/ninjaMesher <caseDir>              # mesh the case
build/ninjaMesher                        # mesh the current directory (like snappyHexMesh)
mpirun --bind-to none -np 4 build/ninjaMesher <caseDir> # parallel run
# (--bind-to none matters: default core binding pins all OpenMP threads
#  of a rank to one core and forfeits the threaded speedup)

build/ninjaMesher --help
build/ninjaMesher --version
build/ninjaMesher --parse-only <caseDir>    # validate the dict, no meshing
build/ninjaMesher --cut-stats <caseDir>     # cut-stage statistics only
build/ninjaMesher --refine-stats <caseDir>  # refinement statistics only
```

## The `ninjaMeshDict`

The dict uses an OpenFOAM-like subset syntax. Minimal annotated example:

```c++
domain
{
    min (-1 -1 -1);      // background box corners
    max ( 1  1  1);
    n   (40 40 40);      // base grid cell counts per axis
}

patches                  // one entry per box side: xmin xmax ymin ymax zmin zmax
{
    xmin { type patch; name inlet;  }
    xmax { type patch; name outlet; }
    ymin { type wall;  name bed;    }
    ymax { type patch; name top;    }
    zmin { type wall;  name side0;  }
    zmax { type wall;  name side1;  }
}

geometry                 // STL surfaces to cut against (paths relative to caseDir)
{
    body.stl { name body; }   // becomes the wall patch name
}

refinementGeometry       // optional: STLs that only SIZE the grid -- never cut,
{                        // never a patch, never layered
    fin.stl { }          // e.g. the pieces of a wall to refine at different levels
}

locationInMesh (0.5 0.5 0.5);   // a point inside the fluid region

refinement               // zero or more named rules
{
    // distance/radius/min/max are ABSOLUTE lengths (mesh units)
    nearBody  { type surface; stl body.stl; distance 0.1; level 2; }
    box1      { type box;    min (0 0 0); max (1 1 1); level 1; }
    ball      { type sphere; centre (0 0 0); radius 0.5; level 1; }
    // `stl` may name a `geometry` or a `refinementGeometry` surface.
    // Several rules may name the SAME stl, which is how you grade bands
    // outward; a cell takes the DEEPEST level of every rule containing
    // it, so levels only ever go up and rule order does not matter.
    bodyFar   { type surface; stl body.stl; distance 0.5; level 1; }
}

layers                   // optional prismatic inflation layers per surface
{
    absorbVolFrac 0.3;   // optional, default 0.3
    body.stl
    {
        nLayers 4;              // default 4
        expansionRatio 1.5;     // default 1.5
        // Thickness of the OUTERMOST layer as a FRACTION of the local
        // cell (snappy's `relativeSizes true` convention). Default 0.3.
        // The total is derived:
        //   total = finalLayerThickness * h * sum_{i<nLayers} ratio^-i
        // where h is the cell size at this surface's refinement level.
        finalLayerThickness 0.4;
        // smoothRadius 1.0;  // landing-surface smoothing radius in LOCAL
                              // GRID-CELL units (x the cell size at the
                              // wall); default 1.0, 0 disables. Keep well
                              // under the smallest feature radius.
    }
}
```

Every layer key has a default, so a bare `body.stl { }` is already a valid
and sensible request (4 layers, ratio 1.5, outermost layer 0.3 of the local
cell).

**Layer thickness is always relative to the local cell.** Absolute
`totalThickness` / `firstLayerThickness` are rejected, with a message that
quotes the equivalent fraction for that entry so an old dict says exactly
what to write. The reason is that layer quality depends on the
thickness-to-cell ratio while the cell size varies with refinement, so one
absolute length silently means different things on different parts of the
same wall — and cannot be right on all of them at once.

**Duplicate keys are a parse error.** Within any one block a key may
appear only once — a repeated region, geometry or layers name is refused
rather than silently overwriting the earlier entry.

**Units.** All coordinates and distances in the dict are absolute (mesh
units): `domain.min/max`, `locationInMesh`, and refinement `distance` /
`radius` / `min` / `max` / `centre`. The grid-relative lengths are
`smoothRadius` (multiples of the local cell size at the wall) and
`finalLayerThickness` (a fraction of it), so both scale with refinement
level.
`expansionRatio`, `absorbVolFrac`, and refinement `level` are dimensionless.

## Diagnostics

By default the mesher prints a banner, a progress line per pipeline stage
(including one per layer-march step), and a final statistics summary.
Useful environment variables (all default off):

- `NINJA_VERBOSE=1` — full statistics tail (landing/thickness counters, …)
- `NINJA_STAGE_TIMES=1` — per-stage wall-clock timing
- `NINJA_MEM_LOG=1` — memory usage logging
- `NINJA_PRELAYERS_DUMP=<dir>` — also write the pre-layer mesh for inspection

## Tests and benchmarks

Both need the OpenFOAM environment sourced (see Dependencies).

```sh
cd build && ctest                    # unit + end-to-end gates (fast)
Benchmarks/run_all.sh                # fast benchmark suite (~1-2 min), needs checkMesh
Benchmarks/run_all.sh --all          # full suite incl. large CAD cases (~1 h)
```

Benchmark gates (checkMesh verdicts, volume/area bands, wall-clock caps) live
in `Benchmarks/gates.toml`; results land in `Benchmarks/results/`.

## How it works

[`DESIGN.md`](DESIGN.md) is the companion document: the pipeline stage by
stage, what each guarantee rests on, and the measured behaviour outside
the comfortable envelope.

## License

Copyright 2026 [Nicolás Diego Badano](https://hydronumerical.com/nicolasbadano/).

NinjaMesher is licensed under the [Apache License 2.0](LICENSE).

