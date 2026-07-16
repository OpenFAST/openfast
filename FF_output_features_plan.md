# FAST.Farm Output Features — Implementation Plan

**Scope.** Two new output features for FAST.Farm's `VISUALIZATION` section,
driven by the "Large Wind Farm Task" (§2.8 and §2.9) and issue
[#2383](https://github.com/openFAST/openfast/issues/2383):

1. **Terrain-following point-cloud sampling** (task §2.8)
2. **Axis-aligned planar sampling** (task §2.9, Snippet 6 syntax)

> Feature 3 (bracket-suffix subregion on the existing `NOutDisWind{XY,YZ,XZ}`
> slices) has been dropped: Feature 2's arbitrary-plane grammar covers the same
> use case (axis-aligned planes with user-controlled extents) without a
> parallel input syntax to maintain.

Build/test harness: `build-docker-single-debug/`. Regression baseline:
[reg_tests/r-test/glue-codes/fast-farm/LESinflow/FAST.Farm.fstf](reg_tests/r-test/glue-codes/fast-farm/LESinflow/FAST.Farm.fstf).
Local test invocation: `ctest -R LESinflow`.

> No code changes yet — this document is the implementation plan. All
> design decisions have been locked in (see §9); implementation can start
> from the commit order in §8.1. Post-implementation activities (docs,
> `build-docker-double` build + regression run, r-test port) are in §10.

---

## 0. Branch context

* **Working branch:** `f/FF_sliceOutput` (fork off `origin/dev` at
  `e4e6e5229` "Update r-test: RM1 floating tank case"). No WIP commits yet.
* **Divergence from upstream `dev`** (as of last fetch — `origin/dev` has
  advanced past this branch by ~40 commits; not yet rebased):
  * Nothing FAST.Farm-slice-related has changed. The only FAST.Farm-adjacent
    diffs are the `NumDFull` / `NumDBuff` default-argument type change from
    `IntKi` → `ReKi` in
    [FAST_Farm_IO.f90:775–783](glue-codes/fast-farm/src/FAST_Farm_IO.f90#L775)
    and a redundant-cast removal in
    [FAST_Farm_Subs.f90:211](glue-codes/fast-farm/src/FAST_Farm_Subs.f90#L211).
    Neither touches the `VISUALIZATION` block, AWAE, or any anchor below.
  * The other ~40 commits are HydroDyn (nonlinear Froude-Krylov), MoorDyn
    Syrope, ServoDyn StC logic, and docs/registry regen — no conflict
    surface for slice output.
  * **Recommendation:** rebase `f/FF_sliceOutput` onto `origin/dev` before
    the "foundation" refactor lands so the shared helpers sit on top of
    current tip. Trivial rebase — expected to be conflict-free.
* **Verified anchors** (still valid on both this branch and `origin/dev`):
  * `maxOutputPlanes = 999` and the three `NOutDisWind{XY,YZ,XZ}` range
    checks live at
    [FAST_Farm_IO.f90:13, 1178–1180](glue-codes/fast-farm/src/FAST_Farm_IO.f90#L13).
  * `VISUALIZATION` block reader begins at
    [FAST_Farm_IO.f90:920](glue-codes/fast-farm/src/FAST_Farm_IO.f90#L920).
  * `ExtractSlice` (axis-aligned trilinear) at
    [AWAE.f90:62](modules/awae/src/AWAE.f90#L62).
  * Slice VTK emit block at
    [AWAE.f90:1992–2050](modules/awae/src/AWAE.f90#L1992).
  * `STRUCTURED_GRID` writer + `.vtk.series` sidecar precedent lives in
    [modules/awae/src/AWAE_vtk.f90](modules/awae/src/AWAE_vtk.f90) —
    Verified 2026-07-16: `AWAE_vtk.f90`, `Write_Planes_Data`, and
    `Write_WakePlane_Series` — which earlier revisions of this plan
    treated as reusable prior art — **are not present** in this
    branch or in `origin/dev`. The XML VTK writers and `.series`
    sidecar helper below are therefore built from scratch, not
    refactored from existing helpers.
* **Branch-name signal:** `sliceOutput` is the natural home for **Feature 2**
  (axis-aligned planar slices). Feature 1 (terrain-following point cloud)
  and the foundation refactor also live on this branch — all three ship
  in a single PR (see §8) rather than split across siblings.

---

## 1. Anchors in the current code

Everything the two features touch lives on the AWAE side. The FAST.Farm glue
code reads the `VISUALIZATION` block into `AWAE_InitInp`, and AWAE writes VTK
at the low-res "skip" rate `p%WrDisSkp1` inside `AWAE_CalcOutput`.

| Concern | File | Key symbol / lines |
|---|---|---|
| Input file `VISUALIZATION` block | [glue-codes/fast-farm/src/FAST_Farm_IO.f90](glue-codes/fast-farm/src/FAST_Farm_IO.f90) | `ReadPrimaryFile` §Visualization (`WrDisWind`, `NOutDisWindXY`, `OutDisWindZ`, …); lines ~920–940 |
| Input validation | [glue-codes/fast-farm/src/FAST_Farm_IO.f90](glue-codes/fast-farm/src/FAST_Farm_IO.f90#L1178) | `ValidateFarmInputData` bounds check `[0, 999]` for `NOutDisWind*` |
| Registry (persistent params + `InputFile` copy) | [modules/awae/src/AWAE_Registry.txt](modules/awae/src/AWAE_Registry.txt) | Lines 36–41 (input file) and companion `ParameterType` block |
| Type file (generated) | [modules/awae/src/AWAE_Types.f90](modules/awae/src/AWAE_Types.f90) | Regenerate from registry |
| Copy input → parameters | [modules/awae/src/AWAE.f90](modules/awae/src/AWAE.f90#L1217) | `AWAE_SetParameters`, `p%NOutDisWindXY = …` |
| Slice validity checks | [modules/awae/src/AWAE.f90](modules/awae/src/AWAE.f90#L1384) | `OutDisWindZvalid` / `OutDisWindXvalid` / `OutDisWindYvalid` |
| Slice allocation | [modules/awae/src/AWAE.f90](modules/awae/src/AWAE.f90#L1473) | `m%OutVizXYPlane`, `m%OutVizYZPlane`, `m%OutVizXZPlane` |
| Slice extraction | [modules/awae/src/AWAE.f90](modules/awae/src/AWAE.f90#L62) | `ExtractSlice` (axis-aligned trilinear along one axis) |
| Slice emit (VTK) | [modules/awae/src/AWAE.f90](modules/awae/src/AWAE.f90#L1990) | `AWAE_CalcOutput` VTK block; uses `WrVTK_SP_header` + `WrVTK_SP_vectors3D` |
| Full 3-D disturbed dump | [modules/awae/src/AWAE.f90](modules/awae/src/AWAE.f90#L1992) | `WriteDisWindFiles` |
| Low-res data source | `m%Vdist_low_full(3, 0:nX-1, 0:nY-1, 0:nZ-1)` (SiKi) | Populated by `LowResGridCalcOutput` |
| Grid metadata | `p%LowRes%oXYZ`, `p%LowRes%dXYZ`, `p%LowRes%nXYZ` | Structured, axis-aligned, uniform |

The existing slice-VTK path (`WrVTK_SP_header` + `WrVTK_SP_vectors3D`) writes
`STRUCTURED_POINTS` — mandatory axis-alignment + uniform spacing. **The new
features cannot use that writer**; they need `STRUCTURED_GRID`
(explicit points, arbitrary orientation) and `POLYDATA` / `UNSTRUCTURED_GRID`
for point clouds. There is a precedent: [modules/awae/src/AWAE_vtk.f90](modules/awae/src/AWAE_vtk.f90)
already writes `STRUCTURED_GRID` for wake planes (`Write_Planes_Data`) and
JSON `.vtk.series` sidecars. We should factor and reuse those helpers.

---

## 2. Feature 1 — Terrain-following sampling (task §2.8, issue #2383 pt. 2)

### 2.1 What the user gets

A cloud of sample points is generated from a **surface** (STL file **or** point
list) plus one or more **offsets along a normal**. At each output step, the
disturbed low-res wind is interpolated at every point and written to a single
per-slice VTK file. The primary use case is a mountain / bathymetry / seabed
surface swept at hub height and rotor tops.

### 2.2 Proposed input syntax

Following pdf §2.8 Snippet 5 with two extensions: (i) list of offsets on one
line (so a single STL becomes multiple sheets), (ii) explicit "cloud file"
alternative to STL for the issue-#2383 use case. Source type is selected by a
`SourceType` flag column (`STL` or `Point`), followed by a single `FileName`
column — cleaner than two mutually-exclusive filename columns.

```
--- TERRAIN-FOLLOWING SAMPLING ---
2                   NumTerrainSlices  - Number of terrain-following sample surfaces (-) [0 to 99]
SliceName    Offsets(m)     OffsetNormal   SourceType   FileName
(-)          (m,list)       (-)            (STL|Point)  (quoted string)
"terr"       100 200 300    default        STL          "terrain.stl"
"hubHt"      140.0          (0 0 1)        Point        "hub_points.txt"
```

* `Offsets` is a *comma- or space-separated* list; each value spawns one
  displaced sheet. Empty / `0` = surface itself.
* `OffsetNormal` — unit vector direction to displace. `default` = STL facet
  normal (per-facet) if `SourceType = STL`; for `SourceType = Point` a
  numeric 3-vector is required.
* `SourceType` — case-insensitive keyword. `STL` reads an ASCII or binary STL
  and uses facet vertices (deduplicated) as sample points, and facet normals
  as the local `+z_local` when `OffsetNormal = default`. `Point` reads a
  plain-text point-cloud file (see below).
* `FileName` — path to the STL or point-cloud file. Relative paths are
  resolved relative to the primary input file (standard FAST.Farm
  convention).

#### 2.2.1 Point-cloud file format (`SourceType = Point`)

Plain-text ASCII, one sample point per line, three floating-point columns
`x y z` in metres in the FAST.Farm global coordinate frame (same frame as
`WT_X`, `WT_Y`, `WT_Z` and the low-res grid origin `X0_Low`, `Y0_Low`,
`Z0_Low`). Parsing rules:

* **Delimiter:** whitespace (any run of spaces or tabs) *or* comma. This
  covers both `.txt` (whitespace) and `.csv` (comma) files with a single
  parser. File extension is informational only — content dictates parsing.
* **Comments:** lines whose first non-blank character is `#` or `!` are
  skipped. Trailing comments (`#`/`!` after the three numeric fields) are
  also skipped.
* **Blank lines:** ignored.
* **Header row:** optional. Any line that fails to parse as three numerics
  is treated as a comment *if* it is among the first non-comment line;
  otherwise it is a hard error. This lets `.csv` exports from GIS tools
  with an `x,y,z` header row load unchanged.
* **Column order:** fixed to `x y z`. No permutation is inferred from
  header names — keeps the parser deterministic.
* **Point count:** derived from the file (no `NumPoints` header). Warn if
  `0` points parsed. No policy cap on `NPts` — the NWTC-library allocator
  reports byte counts and aborts cleanly if the OS refuses the request,
  and no other FAST.Farm sizing input (`NumTurbines`, `NumPlanes`, grid
  dimensions) is capped either. An advisory `ErrID_Info` line is emitted
  at init when a slice's total sample-location count (`NPts × N_offsets`)
  is unusually large — currently ≥ 10 000 000 — stating the count and
  the estimated persistent RAM. See "Memory footprint" below for the
  arithmetic.
* **Duplicate points:** not de-duplicated (unlike STL vertices). The
  user's file is trusted as the intended sample set.
* **No connectivity.** The output VTK is a `POLYDATA` "vertex" mesh — one
  VTK vertex cell per input point.

Example minimal file (both styles are valid input):

```text
# hub-height sample cloud, WGS84 → farm-local (m)
    0.0     0.0   140.0
  250.0     0.0   140.0
  500.0     0.0   140.0
  250.0   250.0   140.0
```

```csv
x,y,z
0.0,0.0,140.0
250.0,0.0,140.0
500.0,0.0,140.0
250.0,250.0,140.0
```

Both files, loaded via `SourceType = Point` with a single `Offsets = 0` and
`OffsetNormal = (0 0 1)`, produce the same four-point VTK output at each
`WrDisDT` step.

#### 2.2.2 Memory footprint and the size advisory

Persistent RAM for one terrain slice, given `N_pts` sample points read from
STL or the point-cloud file and `N_off` offsets on the same line, is
dominated by three arrays sized `(3, N_pts, N_off)`:

| Array | Type | Bytes / location | Purpose |
|---|---|---|---|
| `Pts` (displaced coords) | `ReKi` (4 B) × 3 | 12 | Precomputed once at init, reused every step |
| `Vel` (sampled velocity) | `SiKi` (4 B) × 3 | 12 | Per-step scratch; can be released between steps if memory-tight |
| `Valid` (in-domain mask) | `LOGICAL` (usually 4 B) | 4 | Precomputed once at init |

Total ≈ **28 bytes per sample location** (all three resident during a step).
Plus the source data itself: ~50 B/vertex for an ASCII STL, ~30 B/point for
the plain-text point cloud (read-once, released after parse).

Realistic terrain slices are modest:

| Scenario | N_pts | N_off | Total locations | RAM |
|---|---|---|---|---|
| Seagreen-scale STL at 100 m spacing over 27 km × 20 km | ~ 54 K | 1 | 54 K | ~ 1.5 MB |
| Same, 3 hub-height offsets | ~ 54 K | 3 | 162 K | ~ 4.5 MB |
| Fine mountain STL at 10 m spacing over 20 km × 20 km | ~ 4 M | 1 | 4 M | ~ 112 MB |
| Same, 5 rotor-height offsets | ~ 4 M | 5 | 20 M | ~ 560 MB — triggers advisory |

**No hard cap is imposed.** The NWTC-library allocator (`AllocAry`) already
fails loud with a byte-count error message if the OS refuses the request,
and no other FAST.Farm sizing input (`NumTurbines`, `NumPlanes`, grid
resolutions) is policy-capped either. Adding a `MaxTerrainPoints` knob
would be one more input for users to learn without a corresponding safety
benefit — the OS-level failure mode is already clean.

What we **do** emit is a one-shot `ErrID_Info` line at init when the total
sample-location count for a slice is unusually large (currently
`N_pts × N_off ≥ 10 000 000`, roughly 280 MB persistent RAM). The message
reports the count, the estimated RAM, and the source file — so an accidental
raw-LiDAR-scan-scale input (100 M+ points) surfaces immediately in the log
rather than after a long silent init. The threshold is a compile-time
constant in `AWAE.f90`, easily tuned.

### 2.3 Test geometry

Use [`make_mountain_terrain.py`](https://github.com/user-attachments/files/30064480/make_mountain_terrain.py)
to generate both an STL and a `.txt` point cloud of the same surface, placed
inside `reg_tests/r-test/glue-codes/fast-farm/LESinflow/`. The reg-test
variant exercises **both** `SourceType` paths on the identical geometry.

### 2.4 Algorithm

1. **Init** (once, in `AWAE_Init`):
   * Read the surface; produce a canonical points array
     `Pts(3, 1:NPts)` in the farm-global frame.
   * Optional per-point normal `Nrm(3, 1:NPts)`; else broadcast `OffsetNormal`.
   * For each `Offset(k)` produce `PtsOff(:,i,k) = Pts(:,i) + Offset(k)*Nrm(:,i)`.
   * Mark each displaced point in-domain / out-of-domain relative to
     `p%LowRes%oXYZ` + `p%LowRes%nXYZ*dXYZ`.
   * Emit **one** VTK `POLYDATA` topology file per slice up front
     (mesh doesn't change with time). Time-varying output is per-timestep
     scalar/vector arrays only → small, cheap `.vtu` / `.vtp` field files or
     one `.vtp` per step + a `.vtp.series` sidecar (written by the
     new `Write_VTK_Series` helper introduced in the foundation
     commit).
2. **Per-output-step** (inside the `WriteWindVTK` block of `AWAE_CalcOutput`,
   at the `WrTerrainDT`-derived skip interval):
   * Trilinear-interpolate `m%Vdist_low_full` at each in-domain point.
     Out-of-domain points → IEEE quiet NaN on all three velocity
     components (ParaView masks these automatically). The init step has
     already emitted an `ErrID_Warn` line describing the fraction of
     each slice that falls outside the low-res domain, so users see the
     situation before the run starts.
   * Write `<Root>.TerrSlice.<SliceName>.<Tstr>.vtp` (XML `PolyData`) +
     one appended entry in `<Root>.TerrSlice.<SliceName>.vtp.series`.
     Legacy `.vtk` `POLYDATA` remains available via the section-level
     `VTKFormat` switch.
3. **Finalize**: close series sidecars.

### 2.5 Output format

Preferred: **VTK PolyData (`.vtp`) + `.vtp.series`** — matches ParaView's
native series loader; unstructured / arbitrary point counts natively
supported. Fallback ASCII `.vtk` `POLYDATA`/`UNSTRUCTURED_GRID` for parity
with existing legacy writers.

### 2.6 Pros / cons

| Pros | Cons |
|---|---|
| Single mechanism handles STL, cloud, hub-height sheet, wake-follows-hill | New dependency: STL reader (must implement; no NWTC-lib helper) |
| Cheap runtime — mesh is static, only field values change | Point cloud can be huge (millions of pts). Need memory bound + guardrails |
| VTK output composes cleanly with ParaView's series loader | ParaView users unfamiliar with `.vtp` may expect `.vtk` — support both |
| Reuses trilinear interpolation from `ExtractSlice`-style logic | Requires **new** VTK writer path (existing `WrVTK_SP_*` is uniform-grid only) |

### 2.7 Effort estimate (per pdf: task = $12k)

* STL reader + point-cloud parser: 2 days
* Init-time point processing + validity mask: 1 day
* Per-step trilinear sampler + VTK PolyData writer + series sidecar: 2 days
* Registry / input-file plumbing + docs: 1 day
* Reg-test with `make_mountain_terrain.py`: 1 day

---

## 3. Feature 2 — Axis-aligned planar sampling (task §2.9, Snippet 6)

### 3.1 What the user gets

A list of **axis-aligned** planes (XY, YZ, or XZ) placed anywhere in the
low-res domain, with a user-specified 2-D extent. Grid resolution is the
native low-res spacing (`dX_Low`, `dY_Low`, `dZ_Low`) — same as today's
`NOutDisWind{XY,YZ,XZ}` output — so no `npoints` column is required.
Syntax follows pdf Snippet 6 (OpenFOAM-like), with `center` renamed to
`origin` (the plane's corner, not centre).

Arbitrary (non-axis-aligned) orientations are **not** supported in this
feature. They are covered by a natural future extension: point-cloud
sampling at user-defined locations (see §3.7).

### 3.2 Proposed input syntax (Snippet 6, `origin`-based)

```
--- ARBITRARY PLANAR SAMPLING ---
4                   NumPlaneSlices  - Number of axis-aligned planar slices (-) [0 to 99]
SliceName  origin(m)      normal   extent1(m)  extent2(m)
(-)        (m,m,m)        (-)      (m)         (m)
"T1_0D"    (0    -300 0)  (1 0 0)  600         400        # 600x400 YZ plane at x=0
"T1_1D"    (240  -300 0)  (1 0 0)  600         400        # 600x400 YZ plane at x=240
"hubXY"    (-500 -300 0)  (0 0 1)  2000        800        # 2000x800 XY plane at z=0 ("hubXY" name is user-chosen)
"cross"    (-200    0 0)  (0 1 0)  1000        500        # 1000x500 XZ plane at y=0
```

Semantics:

* `origin` — (x₀, y₀, z₀) of the plane's corner in farm-global coordinates.
  The component along `normal` is the plane location; the other two are the
  origin of the 2-D grid.
* `normal` — unit vector, **must be one of** `(1 0 0)`, `(0 1 0)`,
  `(0 0 1)`. Any other value is a fatal input error (see §3.7).
* `extent1`, `extent2` — positive lengths in metres along the two in-plane
  axes, ordered per the table below. The plane spans `[o_a, o_a+extent1]`
  along axis-1 and `[o_b, o_b+extent2]` along axis-2.
* Grid resolution is **fixed to the low-res spacing** (`p%LowRes%dXYZ`).
  `n1 = nint(extent1/d1) + 1`, `n2 = nint(extent2/d2) + 1`.

| `normal` | Plane | axis-1 (extent1 direction) | axis-2 (extent2 direction) |
|---|---|---|---|
| `(1 0 0)` | YZ | +Y, spacing `dY_Low` | +Z, spacing `dZ_Low` |
| `(0 1 0)` | XZ | +X, spacing `dX_Low` | +Z, spacing `dZ_Low` |
| `(0 0 1)` | XY | +X, spacing `dX_Low` | +Y, spacing `dY_Low` |

The user's `SliceName` is free-form (e.g. `T1_0D`, `hubXY`, `cross`) and
appears in the output filename — not tied to the plane's orientation.

No `offsets` list is provided in this feature; users who want the same
plane repeated at several offsets simply add multiple lines (each with the
appropriate `origin`).

### 3.3 Algorithm

1. **Init** (`AWAE_Init`, once):
   * Parse each slice line. **Reject** any `normal` that is not one of the
     three axis-aligned unit vectors (equality within `1e-6`) with the
     error message quoted in §3.7.
   * Derive `n1`, `n2` from the extents and low-res spacing.
   * Check that the plane fits in the low-res domain. If any corner of the
     requested plane lies outside `[oXYZ, oXYZ + nXYZ*dXYZ]`, emit an
     `ErrID_Warn` line at init that spells out the requested bounds, the
     domain bounds, and the fraction of the plane that will be sampled
     as NaN (see per-step step below). The slice is **not** clipped —
     the full `n1 × n2` grid is emitted so ParaView sees the intended
     footprint; out-of-domain nodes carry NaN velocities. Mark each
     slice as valid / invalid (mirrors the existing `OutDisWindZvalid`
     pattern) so a fully-out-of-domain slice is skipped with a warning.
   * Allocate a per-slice sub-array `m%PlaneSliceBuf(:, 1:n1, 1:n2, k)`
     sized `n1 × n2`. Total per-slice RAM is comparable to today's
     `m%OutVizXYPlane` and typically much smaller (subregion instead of
     full domain).
2. **Per-output-step** (in `AWAE_CalcOutput` `WriteWindVTK` block, at the
   `WrPlaneDT`-derived skip interval):
   * For each slice, call an `ExtractSlice`-family routine with additional
     `i_lo, i_hi, j_lo, j_hi` bounds arguments. Same trilinear
     interpolation along the "thin" (normal) axis as today; then a
     bounded copy into the sub-array — no fancy resampling because the
     in-plane grid *is* the low-res grid.
   * Any node whose `(x,y,z)` falls outside the low-res domain is set to
     `ieee_value(0.0_SiKi, ieee_quiet_nan)` on all three velocity
     components. ParaView masks these automatically.
   * Write `<Root>.Plane.<SliceName>.<Tstr>.vts` (XML `StructuredGrid`) +
     one appended entry in `<Root>.Plane.<SliceName>.vts.series`.
     Legacy `.vtk` `STRUCTURED_POINTS` is also supported via a
     section-level `VTKFormat` switch (default: `xml`).
3. **Finalize**: close each per-slice `.vts.series` sidecar.

### 3.4 Reuse from existing code

Feature 2 uses two shared helpers introduced in the foundation commit:

* **Bounded `ExtractSlice` variant** — today's `ExtractSlice` returns the
  full cross-section; the variant accepts `(i_lo, i_hi, j_lo, j_hi)` and
  populates only the subregion. Every existing call site keeps working
  with default "whole plane" bounds.
* **XML VTK writer + `.vts.series` sidecar** — shared with Feature 1's
  `.vtp` writer through a common `Write_VTK_Series(rootName, ext,
  step_list, time_list, …)` helper (new-from-scratch in `AWAE_IO.f90`;
  no pre-existing series writer in this branch to refactor from). The
  XML file body itself is written by feature-specific routines
  (`WriteVTK_StructuredGrid_2D` for Feature 2, `WriteVTK_PolyData` for
  Feature 1) since the XML schemas differ.

The `STRUCTURED_POINTS` writer already used by the existing slice
output (`WrVTK_SP_*` in NWTC library) is not touched by Feature 2 —
its callers in `AWAE.f90:2004-2050` remain intact and continue to
serve the pre-existing `NOutDisWind*` output when `VTKFormat=legacy`.

### 3.5 Pros / cons

| Pros | Cons |
|---|---|
| Snippet 6 grammar is compact and easy to hand-author | Doesn't cover tilted/rotated planes — those wait for the future point-cloud feature (§3.7) |
| Native low-res resolution matches the existing `NOutDisWind*` output; no surprise interpolation artefacts | Grid is coupled to `dX_Low`/`dY_Low`/`dZ_Low` — users who want a coarser sample must post-process |
| XML `.vts` + `.vts.series` loads directly in ParaView with time-series playback; NaN-masking is native | Adds a new writer (shared with Feature 1) — legacy `.vtk` path remains available |
| Fatal-on-off-axis normal is explicit and forward-compatible — users see the future point-cloud path in the error message | Requires a small parser upgrade for the vector-in-parentheses columns |
| Overlaps in capability with existing `NOutDisWind{XY,YZ,XZ}` — clean deprecation path documented in §4.1 (out of scope here) | |

### 3.6 Effort estimate

(Slightly larger than the last iteration — XML `.vts` output pulls in a
small new writer that's shared with Feature 1.)

* Input parser for the mixed-type multi-column table: 1 day
* Bounded `ExtractSlice` variant + call-site updates: 0.5 day
* Init-time validation, bounds-warning, allocation: 1 day
* Per-step wiring in `AWAE_CalcOutput`, per-section `WrPlaneDT` skip logic,
  NaN out-of-domain fill: 1 day
* XML `.vts` writer + `Write_VTK_Series` helper (shared with Feature 1):
  1.5 days
* Reg-test (`LESinflow_ArbSlice`, three planes + one off-axis fatal):
  1 day
* Registry / input-file plumbing + docs: 1 day

### 3.7 Off-axis normals — explicit fatal + forward pointer

Any `normal` other than `(1 0 0)`, `(0 1 0)`, or `(0 0 1)` (within a
tolerance of `1e-6` per component) is a fatal input error. The exact
error text:

```
FATAL: In slice "<SliceName>", the plane normal (<nx>, <ny>, <nz>) is
not axis-aligned. Only (1 0 0), (0 1 0), and (0 0 1) are supported.
A future feature may allow point-cloud sampling at arbitrary locations.
```

This keeps the input grammar stable when arbitrary-orientation support
lands later — users won't have their old decks silently reinterpret if we
ever widen the tolerance.

---

## 4. Cross-cutting streamlining suggestions

To keep the `VISUALIZATION` section coherent as the two new features land
alongside the existing `NOutDisWind*` slices, propose the following
user-facing consolidation:

### 4.1 New section layout

```
--- VISUALIZATION ---
False   WrDisWind          - Write full disturbed wind (unchanged) (flag)
3.0     WrDisDT            - Time step for the existing NOutDisWind* VTK output (s)
DEFAULT VTKFormat          - "xml" (.vtp/.vts + .series sidecars) or "legacy" (.vtk)  [DEFAULT=xml]

--- VISUALIZATION: AXIS-ALIGNED PLANE SLICES ---   (existing, unchanged for now)
2  NOutDisWindXY  - …
90.0, 87.6   OutDisWindZ  - Z coords
… (YZ, XZ unchanged)

--- VISUALIZATION: AXIS-ALIGNED PLANE SLICES (extent-controlled) ---   (Feature 2, Snippet 6 grammar)
0       NumPlaneSlices  - Number of axis-aligned planar slices (-) [0 to 99]
DEFAULT WrPlaneDT       - Time step for Feature 2 VTK output (s) or DEFAULT [DEFAULT=WrDisDT]

--- VISUALIZATION: TERRAIN-FOLLOWING SLICES ---    (Feature 1, pdf Snippet 5 + point-cloud extension)
0       NumTerrainSlices - Number of terrain-following sample surfaces (-) [0 to 99]
DEFAULT WrTerrainDT      - Time step for Feature 1 VTK output (s) or DEFAULT [DEFAULT=WrDisDT]
```

The two `Num*Slices` counters make it easy to bypass either feature.
Each new feature has its own per-section sampling rate (`WrPlaneDT`,
`WrTerrainDT`) that defaults to the existing `WrDisDT` if the user does
not specify one — useful because point-cloud output can be much larger
than an XY hub-height slice, so users often want to sample it less often.
Internally each is converted to a skip count against the low-res time
step (`p%WrPlaneSkp`, `p%WrTerrainSkp`) mirroring the existing
`p%WrDisSkp1` pattern.

Feature 2 and the existing `NOutDisWind{XY,YZ,XZ}` blocks are
functionally overlapping (both produce axis-aligned uniform-grid VTK
slices) — the new block adds user-controlled 2-D extents that the
existing blocks lack. A future clean-up could deprecate the old blocks
in favour of Feature 2's grammar; that is **not** in scope here (would
break every existing large-farm deck without adding capability).

### 4.2 Shared plumbing

* **XML VTK output for both features (default).** `.vtp` (`PolyData`)
  for Feature 1 point clouds, `.vts` (`StructuredGrid`) for Feature 2
  planes, each with a matching `.series` JSON sidecar. Loads directly
  in ParaView with time-series playback; native NaN masking. Legacy
  `.vtk` remains available via `VTKFormat = legacy` for parity with
  today's outputs.
* **One `.series` sidecar writer.** New
  `Write_VTK_Series(rootName, ext, step_list, time_list, …)` helper
  — built from scratch in `AWAE_IO.f90` (this branch has no existing
  series writer to refactor from).
* **Feature-specific XML file-body writers.** `WriteVTK_PolyData` for
  point clouds (Feature 1) and `WriteVTK_StructuredGrid_2D` for planes
  (Feature 2) — separate routines because the XML schemas differ, but
  they share the same sidecar helper above.
* **One trilinear sampler.** New
  `subroutine SamplePointsFromLowRes(p, V, Pts, mask, out)` used by
  Feature 1. Feature 2 uses the simpler bounded `ExtractSlice` variant
  (still trilinear, but on the aligned grid it collapses to 1-D linear
  interpolation along the normal axis).
* **One STL / point-cloud reader.** Small dependency-free ASCII+binary
  STL reader plus a plain-text `x y z` parser in `AWAE_IO.f90`. Reused
  later for the ground-effect mirror model (task §2.2).
### 4.3 Validation

Existing `[0, 999]` bounds check in `ValidateFarmInputData`
([FAST_Farm_IO.f90:1178](glue-codes/fast-farm/src/FAST_Farm_IO.f90#L1178))
extended to the two new counters (`NumPlaneSlices`, `NumTerrainSlices`).
Each slice line individually checked in `AWAE_Init` and any invalid ones
marked `.false.` in a validity array (mirrors the current
`OutDisWindZvalid` pattern) so partial success is still useful.

---

## 5. Registry & type-file changes

New types (add to [modules/awae/src/AWAE_Registry.txt](modules/awae/src/AWAE_Registry.txt)):

```
# --- Axis-aligned planar slice (Feature 2, Snippet 6) ---
typedef  ^  AWAE_PlaneSliceType   character(64)  Name          -  - -  "..."  -
typedef  ^                        ReKi           Origin       {3} - -  "corner (x0,y0,z0)"  m
typedef  ^                        ReKi           Normal       {3} - -  "axis-aligned unit vector; must equal (1,0,0), (0,1,0), or (0,0,1)"  -
typedef  ^                        ReKi           Extent1       -  - -  "in-plane extent along the first non-normal axis"  m
typedef  ^                        ReKi           Extent2       -  - -  "in-plane extent along the second non-normal axis"  m
typedef  ^                        IntKi          NormalAxis    -  - -  "1=X, 2=Y, 3=Z; derived from Normal at init"  -
typedef  ^                        IntKi          N1            -  - -  "nint(Extent1/d1)+1; derived at init"  -
typedef  ^                        IntKi          N2            -  - -  "nint(Extent2/d2)+1; derived at init"  -
typedef  ^                        LOGICAL        Valid         -  - -  "false if init clipping made the slice degenerate"  -

# --- Terrain-following slice (Feature 1) ---
typedef  ^  AWAE_TerrainSliceType character(64)  Name         -  - -  "..."  -
typedef  ^                        IntKi          SourceType   -  - -  "TerrainSrc_STL=1, TerrainSrc_Point=2"  -
typedef  ^                        character(1024) FileName    -  - -  "..."  -
typedef  ^                        ReKi           Offsets      {:} - -  "..."  m
typedef  ^                        ReKi           OffsetNormal {3} - -  "..."  -
typedef  ^                        ReKi           Pts          {:}{:}{:} - - "..." m
typedef  ^                        LOGICAL        Valid        {:}{:}    - - "..." -
```

Both types embed inside `AWAE_InputFileType` and `AWAE_ParameterType` as
allocatable arrays sized `NumPlaneSlices` / `NumTerrainSlices`.

---

## 6. Testing plan

Base case: [reg_tests/r-test/glue-codes/fast-farm/LESinflow/FAST.Farm.fstf](reg_tests/r-test/glue-codes/fast-farm/LESinflow/FAST.Farm.fstf).
1-turbine LES precursor deck, small & fast, already covered by `ctest -R LESinflow`.

Create two sibling decks under `reg_tests/r-test/glue-codes/fast-farm/`:

| Deck | Feature exercised |
|---|---|
| `LESinflow_TerrainSlice/` | STL + point cloud (Feature 1), regenerated with `make_mountain_terrain.py` at deck setup time |
| `LESinflow_ArbSlice/` | 3 axis-aligned planes with distinct 2-D extents — one XY, one YZ, one XZ (Feature 2). Plus one intentionally off-axis line to exercise the fatal error path from §3.7. |

Each deck reuses the existing `Inflow/` and `WAT_MannBoxDB/` directories from
the base LESinflow test to avoid duplicating turbulence data. Reg-test
validation: byte-compare VTK output (small tolerance for XML timestamps) and
the primary `.out` (should be identical to base — only visualization outputs
differ).

Build/test loop:

```bash
cd build-docker-single-debug
cmake --build . --target FAST.Farm
ctest -R LESinflow                    # baseline
ctest -R LESinflow_ArbSlice
ctest -R LESinflow_TerrainSlice
```

Unit tests (Fortran `pFUnit`): add `awae_test` cases for
* trilinear interpolation of a known analytic velocity field,
* STL reader on the mountain STL (round-trip point count).

---

## 7. Risks and mitigations

| Risk | Mitigation |
|---|---|
| VTK explosion — per-step, per-slice, per-offset files → millions of files on long large-farm runs | Default to XML `.vtp`/`.vts` with **one `.series` sidecar per slice**, not per offset; each new feature has its own `WrPlaneDT` / `WrTerrainDT` skip rate independent of the existing `WrDisDT` so users can sample large point clouds less often than plane slices |
| Memory blow-up on huge point clouds | Rely on `AllocAry`'s existing byte-count-annotated OOM error path (consistent with the rest of FAST.Farm); emit an `ErrID_Info` size advisory at init above ~10 M locations so obvious mistakes surface immediately — no hard cap, no new input knob |
| Line-length ceiling in the input parser | Bump `MaxLineLen` in `FAST_Farm_IO.f90`; already needed for many-turbine decks |
| Registry regen churn (many types + big allocatable ranks) | Do it once, generate `AWAE_Types.f90` in the same commit; verify with `openfast-registry` locally before commit |
| Interaction with domain-decomposition (task §2.1.1) | Design the sampler as a callback that takes `Vdist_low_full` + grid metadata; when partitioned domains land, the callback becomes a gather step. No code duplication needed today. |

---

## 8. Implementation sequencing

**Branch scope (locked in):** everything on the current
`f/FF_sliceOutput` branch, shipping as a **single PR**. No parallel
foundation / terrain-sample branches. The commits below are ordered for
reviewability inside the one branch.

### 8.1 Commit / merge order

1. **Foundation commit** — bounded `ExtractSlice` variant, input parser
   `MaxLineLen` bump if needed, new `Write_VTK_Series` sidecar helper
   in `AWAE_IO.f90` (built from scratch — no pre-existing series
   writer in this branch), and
   the two feature-specific XML body writers
   (`WriteVTK_PolyData`, `WriteVTK_StructuredGrid_2D`). No new user
   knobs, no output changes visible to existing decks. Existing
   `NOutDisWind*` output stays on legacy `.vtk` unless the user opts
   into `VTKFormat = xml`.
2. **Feature 2 (axis-aligned planar sampling, Snippet 6 grammar)** —
   the new `NumPlaneSlices` block, `WrPlaneDT` section-level rate,
   `AWAE_PlaneSliceType` registry entry, off-axis-normal fatal path,
   NaN out-of-domain fill, `LESinflow_ArbSlice` reg-test. Feature 2
   ships first because its grammar and writer path validate the
   foundation work; the trilinear sampler and STL reader stay behind
   until Feature 1.
3. **Feature 1 (terrain-following point cloud)** — STL + point-cloud
   readers, `NumTerrainSlices` block, `WrTerrainDT` section-level rate,
   `AWAE_TerrainSliceType` registry entry, `SamplePointsFromLowRes`
   trilinear helper, NaN out-of-domain fill with init bounds warning,
   `LESinflow_TerrainSlice` reg-test.

Each feature ships behind its own `Num*Slices = 0` default, so partial
adoption of the resulting release is trivial.

---

## 9. Design decisions (locked in)

The questions in the previous revision of this plan have been answered
as follows; recording the outcomes here so they don't get re-litigated
later:

1. **VTK format.** *Answer: add XML `.vtp` / `.vts` + `.series`
   sidecars for both features.* Default `VTKFormat = xml`. Legacy
   `.vtk` remains available (`VTKFormat = legacy`) for parity with
   today's existing `NOutDisWind*` output. Feature 1 uses `.vtp`
   (`PolyData`) with `.vtp.series`; Feature 2 uses `.vts`
   (`StructuredGrid`) with `.vts.series`.
2. **Feature 1 input.** *Answer: yes, support both STL and point-cloud
   files.* Selected via the `SourceType` column (`STL` / `Point`) and
   `FileName`. The mountain-terrain reg-test deck exercises both paths
   on identical geometry.
3. **Feature 1 offsets.** *Answer: yes, multiple offsets per line.*
   `Offsets` is a comma- or space-separated list; each value spawns
   one displaced sheet sharing the source geometry.
4. **Out-of-domain sample handling.** *Answer: IEEE quiet NaN on all
   three velocity components.* At init, an `ErrID_Warn` line is emitted
   for each slice that spills outside the low-res domain, quoting the
   requested bounds, the domain bounds, and the fraction of the slice
   that will be NaN-filled. ParaView masks NaN automatically; users
   see the situation before the run starts rather than being surprised
   by silent zeros.
5. **High-res sampling.** *Answer: low-res only for this task.*
   High-resolution grids are per-turbine and would require a different
   sampler; deferred to a follow-up feature if demand materialises.
6. **Output rate.** *Answer: each new feature gets its own sampling
   rate.* `WrPlaneDT` for Feature 2, `WrTerrainDT` for Feature 1, both
   defaulting to the existing `WrDisDT`. Rationale: point-cloud output
   can be much larger than an XY hub-height slice, so users often want
   to sample it less often. Internally each is converted to a skip
   count mirroring the existing `p%WrDisSkp1`.
7. **Sequencing.** *Answer: Feature 2 first, Feature 1 second.*
   Reflected in the commit order in §8.1.
8. **Branch split.** *Answer: single branch, single PR.* All work on
   `f/FF_sliceOutput`; no `f/FF_slice_foundation` or `f/FF_terrainSample`
   split. The foundation, Feature 2, and Feature 1 are three distinct
   commits on the same branch for reviewability.

---

## 10. Post-implementation activities

Once all three feature commits are in place on `f/FF_sliceOutput`, close
out the branch by walking the three activities below **in this order**.
Each ties directly to a user-provided instruction; the code activity
completes first because it produces the artefacts the other two rely on.

### 10.1 Documentation updates

Touch the four Sphinx sources under
[docs/source/user/fast.farm/](docs/source/user/fast.farm/) that describe
the input file and outputs:

| File | Changes |
|---|---|
| `InputFiles.rst` | New sub-sections after the existing `VISUALIZATION` block (~line 831) documenting the two new `NumPlaneSlices` and `NumTerrainSlices` blocks, `WrPlaneDT` / `WrTerrainDT`, the `VTKFormat` switch, and the `SourceType` / `FileName` / `Offsets` / `OffsetNormal` columns. Include one worked example per feature. |
| `OutputFiles.rst` | New paragraph describing the `.vts` / `.vtp` output naming (`<Root>.Plane.<SliceName>.<T>.vts`, `<Root>.TerrSlice.<SliceName>.<T>.vtp`), the `.series` sidecars, and the ParaView loading workflow. Note that NaN-filled nodes are auto-masked. |
| `ModelGuidance.rst` | Add guidance on when to use Feature 1 vs Feature 2 vs the existing `NOutDisWind*` blocks. Emphasize that Feature 2 is the recommended path for new decks needing user-controlled extents. |
| `FutureWork.rst` | Add "Arbitrary-orientation planar sampling" (the tilted-plane extension that Feature 2's off-axis fatal error message forward-references) as a listed future feature; link to issue #2383. |

Also update the top-level `docs/source/user/api_change.rst` with a
compact "FAST.Farm: new VISUALIZATION sub-sections" entry pointing users
at `InputFiles.rst`.

### 10.2 Build & regression run in `build-docker-double`

`build-docker-double/` does not exist on this branch yet — create it
alongside the existing `build-docker-single-debug/` used during
development.

```bash
cd /software-development/GitHub/openfast-5
cmake -S . -B build-docker-double \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DDOUBLE_PRECISION=ON \
      -DBUILD_TESTING=ON \
      -DBUILD_FASTFARM=ON
cmake --build build-docker-double --target FAST.Farm regression_tests -j
```

The `regression_tests` meta-target pulls in every deck registered via
`ff_regression(...)` in [reg_tests/CTestList.cmake](reg_tests/CTestList.cmake#L440-L446).

Run **only** the FAST.Farm suite (there's no need to spend cycles on the
OpenFAST or module tests for this branch):

```bash
cd build-docker-double
ctest -L fastfarm --output-on-failure
```

> **Naming note:** the ctest label is lowercase `fastfarm`, defined in
> [reg_tests/CTestList.cmake](reg_tests/CTestList.cmake#L442). The user
> request said `ctest -L FAST.farm`, but ctest labels are
> case-sensitive; `FAST.farm` matches zero tests on this repo. Using
> `-L fastfarm` matches all seven currently-registered FAST.Farm cases
> plus the two new ones added in \u00a710.3.

Expected outcome for the existing seven decks:

* `AMReX`, `TSinflow`, `LESinflow`, `TSinflow_curl`, `ModAmb_3`,
  `TSinflowADskSED`, `MD_Shared` \u2014 all pass, byte-identical to the
  gold baselines. The features here ship behind `NumPlaneSlices = 0`
  and `NumTerrainSlices = 0` defaults, so existing decks must produce
  unchanged `.out` / `.outb` output.

If any pre-existing test regresses, the root cause is almost certainly
the bounded `ExtractSlice` variant or the parser `MaxLineLen` bump from
the foundation commit \u2014 both should be exercised without behaviour
change, and a regression there indicates a bug in the shared refactor,
not in the features.

### 10.3 Port new test decks to `reg_tests/r-test/`

During development the two new decks (`LESinflow_ArbSlice/`,
`LESinflow_TerrainSlice/`) live wherever the developer put them.
Porting them into the regression suite is a three-step process because
`reg_tests/r-test/` is a git submodule:

1. **Copy the deck directories into the submodule:**

   ```bash
   cp -r LESinflow_ArbSlice     reg_tests/r-test/glue-codes/fast-farm/
   cp -r LESinflow_TerrainSlice reg_tests/r-test/glue-codes/fast-farm/
   ```

   Each deck reuses `../LESinflow/Inflow/` and `../WAT_MannBoxDB/` via
   relative paths in the `.fstf` (the existing `LESinflow` deck already
   demonstrates this pattern) \u2014 no duplication of turbulence data.

2. **Generate gold baseline files.** Run each new deck once in the
   `build-docker-double` binary and copy the resulting `.out` / `.outb`
   / VTK outputs back into the r-test deck as the "known-good" baseline:

   ```bash
   cd build-docker-double/glue-codes/fast-farm/LESinflow_ArbSlice
   ../../../../glue-codes/fast-farm/FAST.Farm FAST.Farm.fstf
   cp FAST.Farm.out FAST.Farm.T1.outb \
      ../../../../reg_tests/r-test/glue-codes/fast-farm/LESinflow_ArbSlice/
   # repeat for LESinflow_TerrainSlice
   ```

   For the VTK outputs, spot-check with ParaView before adopting them
   as baselines. Byte-comparison of XML `.vts`/`.vtp` files has a
   small timestamp-line tolerance handled by the `openfast-python`
   test harness.

3. **Register in `CTestList.cmake`.** Add two lines after the existing
   FAST.Farm block ([reg_tests/CTestList.cmake:440-446](reg_tests/CTestList.cmake#L440-L446)):

   ```cmake
     ff_regression("LESinflow_ArbSlice"     ""  "fastfarm")
     ff_regression("LESinflow_TerrainSlice" ""  "fastfarm")
   ```

   No extra flags needed \u2014 both decks follow the standard FAST.Farm
   template used by `LESinflow`.

4. **Re-run the FAST.Farm suite** with the new registrations:

   ```bash
   cd build-docker-double
   cmake --build . --target regression_tests
   ctest -L fastfarm --output-on-failure
   ```

   Now nine tests should be listed and all pass.

5. **Commit the submodule bump.** `reg_tests/r-test/` is a submodule;
   the two new deck directories are committed inside it, then the
   outer repository's pointer is updated:

   ```bash
   cd reg_tests/r-test
   git add glue-codes/fast-farm/LESinflow_ArbSlice \
           glue-codes/fast-farm/LESinflow_TerrainSlice
   git commit -m "Add FAST.Farm slice-output & terrain-sampling regression decks"
   cd ../..
   git add reg_tests/r-test reg_tests/CTestList.cmake
   git commit -m "Bump r-test submodule; register new FAST.Farm reg tests"
   ```

   The r-test submodule commit is what carries the deck files; the
   outer commit only advances the submodule pointer and registers the
   tests in `CTestList.cmake`.

### 10.4 Definition of done

The branch is ready for PR review when **all** of the following hold:

* `ctest -L fastfarm --output-on-failure` from
  `build-docker-double/` reports **9/9 passing** (7 pre-existing +
  2 new).
* `docs/source/user/fast.farm/InputFiles.rst`, `OutputFiles.rst`,
  `ModelGuidance.rst`, `FutureWork.rst`, and `api_change.rst` all
  render without warnings under Sphinx (`make -C docs html`).
* Both new r-test decks contain a gold `.out`/`.outb` baseline that
  matches a fresh run byte-for-byte.
* `reg_tests/r-test/` submodule pointer in the outer commit refers to
  a submodule commit that contains the two new deck directories.
* The `f/FF_sliceOutput` branch is rebased on `origin/dev` and has
  three logical commits (foundation, Feature 2, Feature 1).
