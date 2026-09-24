# How NinjaMesher works

A companion to the README. The README tells you how to *run* the mesher;
this document explains the design it implements, why each stage is shaped
the way it is, what it guarantees, and — at equal length — what it does
not.

NinjaMesher is a hex-dominant Cartesian cut-cell mesher that writes
OpenFOAM `polyMesh` directly. It is about 14 000 lines of C++17 and has
no third-party dependency beyond MPI and (optionally) OpenMP.

---

## 1. The thesis

Meshing tools generally treat the CAD surface as a constraint to be
satisfied and cell quality as something to be recovered afterwards. Snap
the cells onto the surface; then smooth, relax and repair until
`checkMesh` stops complaining.

NinjaMesher inverts that:

> **Cell quality is the constraint. Geometric fidelity is the residual —
> the quantity that absorbs the compromise when the two collide.**

The argument for the inversion is that the two error sources are not
comparable in kind:

- A **wall in the wrong place** produces an error that is *bounded by
  the local cell size*, *visible* in the output, and *reducible* by a
  lever the user already has — local refinement.
- A **badly shaped cell** produces an error in every gradient and flux
  the discretisation computes at that cell. It is not bounded by
  anything geometric, it is not visible by looking at the mesh, and no
  amount of refinement elsewhere removes it.

In hydraulic engineering — rivers, spillways, fish passages, bottom
outlets — the "true" wall position is already uncertain to centimetres:
survey error, as-built deviation, sediment, biofouling. Spending mesh
quality to chase millimetres against that background is a bad trade.
That is the domain this design was built for, and the domain in which
the trade-off is being claimed. In aerospace or turbomachinery, where
the surface *is* the answer and it is known to microns, the trade lands
differently and this is the wrong tool.

Everything below is downstream of that one sentence.

---

## 2. The pipeline, inverted

The stages, in order. Stages 1–3 are MPI-parallel over a block
decomposition and their results are gathered to rank 0; stages 4–8 run
serially there, with OpenMP threading the geometric queries.

```
  1. Background grid          perfect axis-aligned hexes, no geometry involved
  2. Octree refinement        2:1-graded, box / sphere / surface-distance rules
  3. Classify + intercept     one in/out verdict per grid vertex,
                              one intercept per cut grid edge
        ── MPI gather ──
  4. Cut                      one cut face per cell, convex cells only
  5. Repair                   sliver merge → coplanar-flap repair → planarize
  6. Layer march              march the cut front DOWN to the true wall
  7. Connectivity drop        discard fluid unreachable from locationInMesh
  8. Write                    polyMesh + cellLevel/pointLevel + wall-distance VTK
```

The shape worth noticing is where the geometry enters. A snapping
mesher starts from the surface and works outward. NinjaMesher starts
from a perfect lattice and works *inward*, and the wall is the last
thing it touches — stage 6, after every quality decision has already
been made.

### The offset trick

For a case with inflation layers, stage 4 does not cut against the STL.
It cuts against the **implicit offset solid**

```
    { p : inside_STL(p) }  ∪  { p : dist(p, STL_s) ≤ t_s  for some surface s }
```

where `t_s` is surface *s*'s total requested layer thickness. The cut
therefore lands one full layer stack *away* from the wall, on a surface
that is smoother than the STL by construction (concave features merge
instead of self-intersecting; convex edges round off at radius `t`).

The cut face produced there is not a wall face — it is the **layer top**,
an internal face. Stage 6 then marches that front down the
distance-to-surface field toward the true wall in graded steps, emitting
one prism per front face per step.

So the boundary layer is not inflated *out of* a snapped skin, with the
core mesh pushed aside to make room and quality crossed off afterwards.
It is built *into* a void that was reserved for it before any cell was
cut, out of a front that was already clean.

Crucially there is **no projection map** — no "find each mesh point's
nearest wall point and pull it there". That map is discontinuous across
concave features and across the medial axis, and repairing those
discontinuities is precisely the tune-the-mesh-against-exact-geometry
loop this design exists to avoid. Geometry enters the march only as a
scalar field the front flows through.

---

## 3. Stage by stage

### 3.1 Grid and refinement

A background box of `n = (nx, ny, nz)` hexes, then octree refinement
driven by named rules: `box`, `sphere`, and `surface` (refine within a
distance *d* of a named STL). A cell takes the deepest level of every
rule that contains it, so levels only ever go up and rule order is
irrelevant — several bands over the same STL is how you grade outward.

Two invariants are enforced here and relied on downstream:

- **2:1 grading.** Adjacent cells differ by at most one level. This
  bounds the volume jump across every internal face, and it is exactly
  the condition OpenFOAM's AMR (`dynamicRefineFvMesh`) requires — which
  is why `cellLevel` and `pointLevel` are written alongside the mesh
  and a NinjaMesher mesh can start an AMR run directly.
- **Refine, then cut.** The interface polyhedra at level transitions
  exist *before* the cut, so the cutter never sees a special case: it
  operates on general convex polyhedra, of which a hex is the 6-face
  instance.

The surface-distance rule is inherently iterative — splitting a cell
exposes children whose distances only then become meaningful — so the
mark/split/re-grade loop runs to a fixpoint, with the 2:1 halo exchanged
across MPI ranks until an `Allreduce` says no rank moved.

### 3.2 The cut: shared intercepts and parity ("A+parity")

This is the load-bearing idea, and it is a topology decision rather than
a geometry one.

Every quantity that two neighbouring cells could disagree about is
computed **exactly once, at the entity that owns it**, and then read by
every cell that touches it:

- **Vertex in/out** — one ray-parity test per *grid vertex*, referenced
  against `locationInMesh` so the STL's winding convention never
  matters.
- **Edge intercept** — one intersection per *grid edge*, stored and
  reused bit-identically by all four (or more) cells on that edge.

A cell's cut is then a pure function of that shared data: keep the
fluid-side vertices, keep the kept-side stubs of cut edges, close the
cell with a single face through the shared intercepts.

**Watertightness is therefore true by construction, not defended with
tolerances.** Two cells cannot disagree about a shared vertex or a
shared intercept because there is only one of each. Compare the
alternative — fitting a plane per cell and clipping against it — where
neighbouring cells fit *different* planes and a near-surface vertex can
be classified two ways. That is the classic silent watertightness
killer, and it is structurally impossible here.

A best-fit plane (AMReX-style, from the divergence of face-area
fractions) does survive, but only as the cut face's *orientation hint*.
It never classifies anything.

The price, stated plainly: **one cut face per cell**. A thin plate or a
narrow gap smaller than a cell needs two cuts in that cell and gets one
— or, if both its faces cross every grid edge an even number of times,
gets *zero*, and the feature is simply not in the mesh — which is
what asking for cells larger than the feature means. There is a
benchmark, `bm_thin_plate`, whose entire purpose is to demonstrate this
failing loudly enough that nobody discovers it in production. The
user's lever is refinement, never a more elaborate cutter.

### 3.3 The half-grid

The offset cut cannot be taken against the grid alone. When the
requested thickness `t` is a multiple of the cell size `h`, the level
set `d = t` lands *exactly on grid planes*, and `dist ≤ t` is decided by
floating-point noise. Two vertices that are geometrically mirror images
land on opposite sides of the comparison:

```
  |0.3 − 0.25| = 0.04999999999999999   ≤ 0.05  → TRUE
  |0.75 − 0.7| = 0.050000000000000044  ≤ 0.05  → FALSE
```

Two ULPs apart, in opposite directions, because 0.3 and 0.7 are not
binary-exact. An offset region decided that way runs one grid plane
deeper on one side of a symmetric structure than the other, and the two
halves get structurally different layer fronts.

An epsilon does not fix that. Instead the surface is snapped to the
**nearest node of the half-grid** (grid vertices ∪ edge midpoints), and
each vertex is classified into three states against `φ(v) = d(v) − t`:

```
   |φ| < band  →  ON        band = h_local / 4
    φ < 0      →  IN
    φ > 0      →  OUT
```

`band` is not a tolerance: `h/4` is the Voronoi radius of a half-grid
node — the exact width at which "which node is nearest" changes answer.
Vertices whose φ is pure noise (φ = −1.4e−17 and +4.4e−17 on the
rotated-cube case) land inside the band and are both ON, whichever side
of zero the noise falls on. Symmetry is preserved by construction rather
than by a tie-break that a compiler flag could flip.

`h_local` is per vertex: the size of the finest cell touching it, read
off the graded octree, so the verdict stays a pure function of position
(one answer per vertex on every rank). One band sized from the finest
level in the whole case is not equivalent: on a coarser wall it sends
almost every crossing to the edge midpoint, up to `h/2` off, and the
layer stacks there come out too thin to keep (measured: 19474 → 10601
dropped wall faces on `bm_layers_wfp` from this alone).

Quantisation error halves as a side effect: reachable cut positions go
from spacing `h` to `h/2`, so worst-case wall displacement goes from
`h/2` to `h/4`.

### 3.4 Repair, in three passes

The cut leaves three specific, understood defects, each with its own
pass. All three are **refuse-rather-than-damage**: when the repair
cannot be made validly, the pass declines and leaves the cell alone.

1. **Sliver merge.** The offset cut makes its slivers exactly at the
   layer top — the worst possible place, since that is where velocities
   are highest and a tiny cell there pins the solver time step. A
   sub-threshold cut cell (`absorbVolFrac`, default 0.3 of the cell) is
   merged into a kept neighbour rather than removed. A merge is refused
   if the result would be non-manifold, non-convex (centroid outside one
   of its own facets — checkMesh's face-pyramid criterion), fail the
   transcribed per-face geometry tests, or exceed a face-skewness bound
   of 3.0 — on its internal faces and, with checkMesh's boundary
   formula, on its wall faces. Each refusal reason is counted and
   printed.

2. **Coplanar-flap repair.** A chord chain that grazes a foreign ON
   sheet can leave a face parented to the wrong side. This is a
   topology repair with an exact closure identity, not a tolerance, so
   it is always on.

3. **Planarize.** A warped layer-top face is split into a fan. The
   tolerance is derived, not configured: `0.1 × t₀`, a tenth of the
   thinnest wall layer — the length scale the layer top actually has to
   resolve. The warp distribution is bimodal (numerical noise near zero,
   genuine warp far above), so every tolerance inside the gap gives a
   bit-identical mesh; rules based on `dx` or `h_local` land past the
   gap edge and are measurably wrong. A fan is accepted only if the
   owner cell, re-centred by the split, still passes the per-face
   geometry tests and every fan triangle stays under checkMesh's
   boundary skewness; otherwise the face ships warped but valid.

### 3.5 The march

For each layer step, every front point:

- **aims** along `normalize(closestPoint(p) − p)`, the descent direction
  `−∇d`, computed from the closest point rather than a finite-difference
  gradient (which is noisy on a faceted STL), then smoothed over the
  front's point graph;
- **steps** a graded fraction of *its own remaining distance*, so
  grading is expressed in field units and a point clamped on an earlier
  step automatically re-aims with whatever distance it has left;
- **is clamped** by minimum resulting layer height, maximum angle
  against the local front normal, maximum tangential stretch, prism
  non-inversion, and — transcribed from OpenFOAM's own
  `primitiveMeshTools` — checkMesh's **aspect-ratio** metric against its
  own default ceiling of 1000. A prism that survives all of that is still
  refused if one of its faces is **severely non-orthogonal** (checkMesh's
  70°): where an oblique wall meets a grid-aligned one the wall faces are
  needles, and prisms grown from needles meet on their short sides with
  their centres displaced along the strip.

The march runs **outermost layer first**; the last step lands on the
wall. Thicknesses are geometric, thinnest at the wall:
`t_j = t₀·r^j`, `t₀ = t / Σ r^j`.

Two cross-step behaviours are deliberate and are what "graceful
degradation" actually means here:

- A face whose prism fails validation at step *k* **reverts to a plain
  wall face** and rejoins the wall bucket, so step *k+1* re-marches it
  with the distance it has left. The result is locally *fewer, thicker*
  layers — not a hole.
  The exception is a face that shares a front point with a stack that
  *did* step: that point has already been extruded, and a front point
  is extruded at most once, so the face is **held** as a wall face for
  the rest of the march. The same rule holds the terrace's seam faces
  (the side quads of the kept prism). Without it the shared point would
  march twice down the same line and land twice on the same wall spot —
  duplicate points and a zero-width crack.
- A point frozen by the non-inversion clamp re-aims on the next step for
  the same reason.

Where the march declines entirely, the wall simply stays at the offset
cut — away from the true wall, by a distance the mesher reports
(`infeasibleWallArea`). The mesh is valid and honest about being coarse
there.

**Landing smoothing.** With `smoothRadius > 0` (default 1.0, in local
cell sizes) the front lands on a mollified version of the distance
field. This is what makes near-complete marches possible on real CAD —
but it also means sharp creases round at the smoothing-radius scale, and
features comparable to that radius can shrink or seal. Everything
smoothing changes is measured: `smoothMovedArea`, and per-region
`smoothSealedRegion*` records where the smoothed and exact fields
disagree on *sign* — i.e. exactly where a passage sealed or a rib
evaporated. A **seal guard** refuses to emit a wall face that landed
inside the solid; the refused area is reported, never quietly written.

The sweep that fixed the default:

| `smoothRadius` | 0 | **1.0** | 1.5 | 2.0 |
|---|---|---|---|---|
| stack completion | 99.00 % | **99.91 %** | 99.91 % | 100 % |
| dropped faces | 1917 | **182** | 182 | 0 |
| sealed wall area | 0 | **0.34 %** | 1.10 % | 1.92 % |

`r = 1.5` buys nothing over 1.0 and seals 3.2× the area; `r = 2.0`
completes the march but seals 1.9 % of the wall — and sealing passages
changes flow topology, which is itself a robustness defect, not a win.

### 3.6 Connectivity, and what gets discarded

A flood fill from `locationInMesh` across internal faces keeps one
connected component. Everything else is fluid the mesh will not
describe, and it is disclosed rather than silently dropped: per-component
cell count, volume and bounding box on stdout, plus a warning on stderr.

This matters more than it sounds. Coarsening a grid until thin plates
stop being resolved does not produce an error — it produces a *valid
mesh of a different flow domain*, where the fill has run through a wall
and merged two compartments. The `cutDiscarded*` lines are how that
shows up.

### 3.7 Parallelism

Refinement and classification are MPI-parallel over a block
decomposition; the cut and layer tail run on rank 0, with OpenMP
threading the geometric queries.

The output is **partition-invariant by construction**, and this is a
design property rather than a tuned one. The parallel stages compute
values that are pure functions of `(position, replicated STL)`, keyed by
exact integer lattice indices — never by tolerance matching — so a
duplicate arriving from two ranks is bit-identical, which the
reconstruction asserts rather than reconciles. The serial tail then runs
on the same reconstituted input regardless of rank count. The
`np_invariance` gate byte-compares every `polyMesh` file at np = 1, 2
and 4 across eight cases; any difference fails, with no tolerance
comparison anywhere.

On measured cost: OpenMP takes the layer tail on the fish-passage case
from 1023 s to 103 s and the whole run from 18:05 to 2:31, bit-identical
either way. MPI cannot reach that stage at any rank count without a
redesign. Recommended invocation is therefore one rank per node with
threads free — and `--bind-to none`, or the default core binding pins
every thread of a rank to one core and forfeits the speedup entirely.

---

## 4. What is guaranteed, and by what mechanism

| Guarantee | Mechanism | Where it is checked |
|---|---|---|
| Watertight topology | Shared per-vertex verdict + per-edge intercept; neighbours cannot disagree | `check_integrity.py` on every case; `expected_fail` never waives it |
| Convex cut cells | Clipping a convex polyhedron with one cut face | Cut-stage post-condition; failures removed, counted, printed |
| Bounded volume jump | 2:1 octree grading | `--refine-stats` `maxAdjacentLevelDiff == 1` |
| No cell worse than checkMesh's own thresholds | checkMesh's per-face geometry tests transcribed from `primitiveMeshTools` and run *as a post-condition* — not paraphrased | Cut post-condition, merge refusal, and in-march prism validation |
| No prism above the aspect-ratio ceiling | checkMesh's `cellClosedness` aspect metric vs its own default of 1000 | In-march validation; failing stack condemned, wall stays at the offset cut |
| No prism face severely non-orthogonal | Angle between a prism face's normal and the line of centres vs checkMesh's own 70° — exact on the face to the cell above, one-sided on side faces whose neighbour is not decided yet | In-march validation. Measured on three windows: 148, 17 and 7 faces above 70° (worst 82.9°) → none, for 0.3 % fewer prisms |
| No wall face inside the solid | Seal guard | Refused area reported per case |
| Partition invariance | Exact integer keys; pure-function values; identical serial tail | `np_invariance`, byte-compare at np = 1/2/4 |

The recurring pattern is worth naming: rather than *paraphrasing* what
checkMesh would object to, the mesher **transcribes checkMesh's own
tests** and runs them as internal post-conditions. A cell that would
fail is never written in the first place.

## 5. What is approximated, and by how much

| Approximation | Bound | Lever |
|---|---|---|
| Wall position | ≤ ¼ of the local cell (half-grid snapping) | Local refinement |
| Sub-cell features | Not represented at all | Refinement; `bm_thin_plate` is the worked example |
| Sharp creases under smoothing | Rounded at the smoothing radius (× local cell) | `smoothRadius 0`, or refine |
| Narrow gaps under the offset cut | Gaps < 2·t seal; no front forms inside | Smaller `t`, or refine |
| Layer depth where infeasible | Locally fewer/thicker layers, or none | Reported as `infeasibleWallArea` |

The sub-cell row is the one to understand before using the tool, and it
is a *consequence of choosing a cell size*, not a defect in the cut. A
feature thinner than a cell is below the resolution the dict asked for.
One cut face per cell cannot represent two surfaces crossing one cell, so
such a feature is not approximated coarsely — it is absent. In the
extreme (`bm_thin_plate`: a plate 0.4 cells thick, centred off any
grid-vertex plane) both its faces cross every grid edge an even number of
times, the cut finds nothing, and the mesh is a valid model of a domain
with no plate in it.

Nothing is reported, because nothing was seen. That follows from the same
choice: the cut asks its questions at grid vertices and grid edges, and a
feature that misses all of them does not exist as far as the mesh is
concerned — exactly as a 1 m grid does not resolve a 1 cm ripple. The
lever is the same one the whole design rests on: refine what you care
about, and the feature appears. Resolution is the user's decision, and
this is what that decision means.

---

## 6. How the claims are verified

### The checkMesh whitelist, stated precisely

The benchmark suite runs `checkMesh -allGeometry`, which is
*stricter* than plain `checkMesh`, and then accepts five specific error
classes:

- `Concave cells (using face planes)`
- `Faces with small interpolation weight`
- `Faces with small volume ratio`
- `concave angles between consecutive edges`
- `Cells with small determinant`

These five are exactly the checks that `-allGeometry` adds and plain
`checkMesh` does not run (`checkConcaveCells`, `checkFaceWeight`,
`checkVolRatio`, `checkFaceAngles`, `checkCellDeterminant` in
OpenFOAM's `checkGeometry.C`). So the honest statement of the result is:

> **Plain `checkMesh` returns `Mesh OK.`** Under the stricter
> `-allGeometry`, what remains is confined to five documented
> quality-floor classes and nothing else.

They are accepted because they are direct consequences of two design
decisions, not latent bugs: the cut face is allowed to be mildly
non-planar, and refinement-interface polyhedra are non-orthogonal by
construction. A case with no STL at all (`bm_diag_jump`) trips the
concave-cells check purely from refinement interfaces, which is the
cleanest demonstration that these are structural rather than cut
defects. Every occurrence is still recorded verbatim in each case's
result file. Any *other* `***` error fails the gate.

### The gates

- **`ctest`** — 38 gates: unit tests plus end-to-end cases asserting
  watertightness, grading, layer grading geometry, np-invariance, RSS
  scaling, dict validation, and a set of real-CAD sub-box windows that
  reproduce specific historical defects.
- **`Benchmarks/run_all.sh`** — 38 cases, 37 of them gated by
  `gates.toml`: checkMesh verdict, analytic volume and wall-area bands,
  ceilings on non-orthogonality / skewness / aspect ratio / wall-clock,
  floors on achieved layer height and minimum cell volume, and caps on
  sealed area and zero-face patches. 14 are marked `slow` and run only
  under `--all`; the rest are the fast set the `benchmarks` ctest gate
  uses.
- **Solver gates** — seven cases run a real solver on the mesh they just
  built: `interFoam` on the four real-CAD cases, plus `simpleFoam` and
  `icoFoam` on synthetic ones. They gate on reaching the requested
  `endTime`, peak Courant number, bounding-line count and cumulative
  continuity error.

That last family exists because **`checkMesh` is not a sufficient
oracle**. It answers questions about geometry and topology — is this cell
closed, is this face skewed, is this volume positive — and it answers
them one cell at a time, against fixed thresholds. It has no notion of a
time step, so it cannot ask the question that decides whether a mesh is
usable: what does this cell do to the *global* solution?

The clearest instance is the small cell. A solver's Courant number goes
as `0.5·dt·Σ|φ|/V`, so a single cell with a small `V` pins the time step
for the entire domain — while checkMesh, comparing that same volume
against its own floor, reports "Cell volumes OK". The cell is valid. The
run is crippled. Nothing in a checkMesh report distinguishes the two
cases, which is why every real-CAD mesh in this suite is also required
to *run*.

---

## 7. Behaviour outside the comfortable envelope

The design claim is not "the mesh is always good". It is:

> Under-resolution is acceptable: the feature stops being represented, a
> thin plate vanishes, a passage closes or leaks, layers get fewer or
> disappear, the wall is approximated more coarsely — the mesh stays
> valid and stays away from the wall. A quality failure is never
> acceptable, at any cell size, for any thickness request.

That was tested by sweeping four real-CAD cases over grid factors ×2,
×4, ×8 and layer-thickness multiples, with layers on and off.

**It holds on all four cases to ×8**, layers on and off — 24 of 24
meshes clear checkMesh's whitelist — at dx 0.96 m on a case designed
at 0.12 m. What degrades is what should: layers retreat from the wall
monotonically and loudly (`infeasibleWallArea` 0.7 % → 17.0 % → 31.3 %
of the wall as the grid coarsens ×2 → ×4 → ×8, counted on the fluid
that is written), the seal guard refuses progressively more landings
(1.5 → 16.8 → 70.6 m²), non-orthogonality stays under 70° at every
factor because the march refuses the prisms that would breach it, and the outer
quality gate stops converging — while the mesh stays valid, because the
in-march guard reverts a failing prism *before* emission. Non-convergence
costs layer coverage, not validity.

**Boundary faces need a check of their own.** checkMesh scores a
*boundary* face against its owner's centre alone, so a well-shaped face
in the wrong place passes every volume and degeneracy test — only
skewness sees it. The shape that exposes this on `bm_layers_wfp` at grid
×4 is a 16 × 20 mm remnant of a grid face, split by hanging nodes, on a
level-2 cell cut down to a 0.9 % wedge: a boundary face only because the
cell across it is removed. Three stages create or re-centre boundary
faces, and each asks checkMesh's own question of them: the cut removes a
cell whose shipped boundary faces breach it (one cell in 25 103 at ×4),
the sliver merge scores a group's wall faces against the merged centre,
and planarize rejects a fan whose triangles breach where the whole
polygon did not.

Doubling layer thickness past the cell size, by contrast, costs almost
nothing: at t/h = 1.93 on the fish passage, max skewness *improves*
(2.646 → 2.605), max aspect ratio is 27.4, and 1.5 % of the wall ends up
without full-depth layers. What breaks prisms is height against *width*,
and it is the grid that moves that ratio — not the thickness on its own.

And the honest counterpart: below the resolution at which a case's
defining feature survives, you get a valid mesh **of a different
geometry**. On the spillway case, fluid volume grows 61 % from ×2 to ×8
as the radial gates' 198 mm skinplates stop being resolved and the flood
fill counts sealed gate bays as open channel. The mesh is fine. It is
just not a model of that spillway any more, and `wallArea`,
`infeasibleWallArea` and the `cutDiscarded*` lines are how it says so.

---

## 8. Where to look in the source

| Concern | File |
|---|---|
| Pipeline orchestration, dict schema, CLI | `src/main.cpp` |
| Grid, octree refinement, 2:1 grading | `src/Refine.cpp` |
| Vertex classification, edge intercepts, half-grid | `src/CutData.cpp` |
| Cut, sliver merge, flap repair, planarize | `src/Cutter.cpp` |
| Layer march, smoothing, prism validation | `src/Layers.cpp` |
| Closest-point / distance queries, AABB bins, BVH | `src/Geometry.cpp` |
| OpenFOAM `polyMesh` output | `src/FoamWriter.cpp` |

The headers carry the design rationale; the `.cpp` files carry the
mechanism. Where a constant has a value rather than a derivation, the
comment next to it says which measurement chose it.
