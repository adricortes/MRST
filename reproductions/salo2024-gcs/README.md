# Reproduction: Saló-Salgado et al. (2024)

This folder reproduces results from: Saló-Salgado, L., Møyner, O., Lie, K.-A.,
Juanes, R. "Three-dimensional simulation of geologic carbon dioxide
sequestration using MRST". Advances in Geo-Energy Research, 14(1): 34-48, 2024.
DOI: 10.46690/ager.2024.10.06.

The input decks come from the `mrst-adblackoil-gcs` repository. Those decks
are equal to the published CO2LabMIT dataset (a copy is in `data/`). The
maintained example code lives in the `co2lab-mit` module of MRST. The scripts
here are headless drivers based on those examples.

## Environment

- GNU Octave 8.4 on Linux (no MATLAB). Graphics toolkit: gnuplot (default;
  see the notebooks section below for when `fltk` is needed and why it
  must stay opt-in per cell, not a session-wide default).
- MRST development version, commit 8aaf7d2.
- Date of the original headless runs: 2026-08-19. Notebooks and the
  well-index appendix added 2026-08-20.

## Contents

- `scripts/` — headless driver and post-processing scripts, plus
  `computeWellIndexPalagiAziz.m` and `find_fault_cell.m` (see
  "Well-index investigation" below).
- `scripts/shim/` — a minimal `delaunayTriangulation` class for Octave.
- `notebooks/` — four Jupyter notebooks reproducing the same results
  interactively (see "Jupyter notebooks" below).
- `figures/` — output figures (PNG).
- `results/` — saved states (`.mat`, Octave v7 format).
- `logs/` — full run logs.
- `data/` — the CO2LabMIT dataset zip from SINTEF, for reference.

## Results

### 1. PVT tables (paper Fig. 2)

Script: `scripts/run_pvt_props.m`. Figure: `figures/fig_pvt_props.png`.
The MRST fluid object reproduces the tabulated output of
`pvtBrineWithCO2BlackOil` to machine precision (max relative difference
about 2e-16) for CO2 density, CO2 viscosity, brine density, and brine
viscosity. Octave cannot run the generator function itself, because that
function needs the MATLAB `table` class. The comparison therefore uses the
shipped output tables.

### 2. PUNQ-S3 relative permeability hysteresis (paper Fig. 4)

Script: `scripts/run_punq_hyst.m` (env vars: KRHYST, NSTEPS, OUTTAG).
Post-processing: `scripts/postproc_punq.m`. Figure:
`figures/fig_punq_sg500y.png`. Runtime: 33 s and 51 s for 100 steps.

| Metric (t = 500 y)      | No hysteresis | Hysteresis |
|-------------------------|---------------|------------|
| Max Sg                  | 0.690         | 0.641      |
| Cells with Sg > 0.05    | 49            | 138        |
| Cells with Sg > 0.3     | 36            | 23         |
| Mean Sg where Sg > 0.05 | 0.505         | 0.192      |

At the end of injection (t = 10 y) both cases are almost identical. This
matches the paper: hysteresis acts during post-injection imbibition, traps
CO2 residually along the migration path, and shrinks the gas cap.

### 3. Molecular diffusion (paper Fig. 6)

Script: `scripts/run_diffusion.m` (env vars: DIFFD, MESHRES, OUTTAG).
Post-processing: `scripts/postproc_diffusion.m`. Figure:
`figures/fig_diffusion_rs30d.png` (coarse mesh, 8,206 cells, 167 s per run).
On the coarse mesh the diffusion effect is small, because numerical diffusion
dominates at 1 cm cells (the paper's Section 4.1 explains this). The
fine-mesh runs (105,279 cells, `figures/fig_diffusion_fine30d.png`,
2,751 s and 3,532 s) resolve the convective fingers: with diffusion the
fingers are wider and farther apart, and the mean dissolved ratio rises
from 2.97 to 3.54 Sm3/Sm3. This reproduces the paper's Fig. 6(a,b).

A full LaTeX report (theory, MRST practice, all figures) is in
`report/report.tex` and `report/report.pdf`.

### 4. Johansen case study (paper Fig. 7)

Not attempted. The model has 90,000 cells over 500 years and needs MEX and
AMGCL acceleration, which Octave does not provide.

## Jupyter notebooks

`notebooks/` has four executable notebooks, one per section of
`report/report.tex`, running MRST through an Octave Jupyter kernel with
markdown narrative paraphrased from the report:

- `01_trapping.ipynb` — structural trapping (Section 3). Cheap; always
  runs live.
- `02_pvt.ipynb` — PVT verification (Section 4). Cheap; always runs live.
- `03_hysteresis.ipynb` — PUNQ-S3 hysteresis (Section 5).
- `04_diffusion.ipynb` — molecular diffusion (Section 6), plus a mesh
  preview (compartment overview and a zoomed-in view with real cell
  edges showing the PEBI polygons) and the well-index appendix described
  below.

The hysteresis and diffusion notebooks default to loading the cached
`results/*.mat` states rather than re-simulating (`RUN_LIVE` /
`RUN_LIVE_COARSE` / `RUN_LIVE_FINE` flags near the top of each, default
`false`); both the cached and live-rerun paths were verified to reproduce
the numbers in this README and the report exactly. Every notebook also
writes its figures back to `figures/`.

Building them surfaced one further portability pitfall, specific to the
interactive kernel rather than to headless scripts: `fltk` is needed for
the 3-D patch plots (`plotGrid` on a 3-D grid) that gnuplot cannot draw,
but setting it as the *session-wide default* graphics toolkit (e.g. in
`~/.octaverc`) silently renders every `imagesc` plot — the type used
throughout this reproduction's own figures — as a blank raster, with no
error. The fix: keep gnuplot as the default and switch to `fltk` only in
the specific cell that draws a 3-D patch plot, immediately followed by
`set(0, 'DefaultFigureCreateFcn', @(src,~) set(src,'visible','on'))` (an
Octave Jupyter kernel forces new figures invisible for inline capture,
which `fltk`/OpenGL cannot print from).

## Well-index investigation

`04_diffusion.ipynb`'s appendix audits the well model MRST uses for the
diffusion case's injector, which sits on a PEBI mesh, not a Cartesian
grid. MRST's `computeWellIndex`/`cellDims` is Peaceman-type but derives
the equivalent radius from the well cell's *axis-aligned node bounding
box*, regardless of its true polygon shape — it never looks at the cell's
actual neighbors or connection transmissibilities. Palagi and Aziz (SPE
22889, 1991; SPE 24072, 1994) derive a well index directly from those
connections instead, which is defined for any polygon and reduces to
Peaceman's formula on a square Cartesian cell.

- `scripts/computeWellIndexPalagiAziz.m` implements Palagi & Aziz's
  simplified model (SPE 24072, Eq. 17), reusing MRST's own `cellDims` and
  `computeTrans` so it matches MRST's kh, skin, and well radius exactly and
  differs only in how the equivalent radius is derived. `addWell` accepts
  a precomputed well index directly, so no MRST internals needed
  modifying.
- `scripts/find_fault_cell.m` ranks the coarse mesh's fault-strip cells
  by bounding-box aspect ratio to find a deliberately elongated test cell.

Results (coarse mesh, well radius 1 mm, skin 0):

| Location | WI Peaceman | WI Palagi-Aziz | Difference |
|---|---|---|---|
| Reservoir well (original, aspect ratio ~1) | 9.034e-14 | 8.157e-14 | -9.71% |
| Elongated fault-strip cell (aspect ratio ~2.4) | 5.555e-16 | 4.583e-16 | -17.50% |

The well-index gap grows with cell irregularity, as expected. Because
this well is rate-controlled, the well index only sets the drawdown
needed to deliver the fixed rate — it never changes how much CO2 enters
the reservoir, so rerunning the coarse-mesh, D=0 case with each well
index changes the dissolved/free-gas fields only at the level of solver
noise, but the bottom-hole pressure gap grows about 230x from the
reservoir well (0.35 Pa) to the fault-strip cell (80.0 Pa) — real, and
tracking the larger well-index gap and the fault's lower permeability
(10 mD vs. 1000 mD), but still negligible against the ~9.9 MPa reservoir
pressure at this problem's flow rate. A pressure-controlled well, a
higher rate, or a lower-permeability target would be expected to show a
larger effect.

## Octave portability notes

1. `modules/upr/util/sortEdges.m` had a malformed block comment. Octave
   parsed the function body as a comment and returned grids unsorted, which
   silently inverted extruded PEBI cells. Fixed in this MRST checkout;
   the fix deserves an upstream pull request.
2. Octave lacks `deval`, `pdist2`, `griddedInterpolant`, and
   `delaunayTriangulation`. The drivers replace the first two inline, put
   `core/utils/octave_only/mrst` ahead on the path for the third, and use
   `scripts/shim/` for the fourth.
3. `matlab_bgl` MEX binaries do not load in Octave. This blocks only the
   cell-based `trapAnalysis` method, which these examples do not use.
