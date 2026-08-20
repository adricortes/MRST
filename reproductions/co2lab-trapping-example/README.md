# co2lab trapping example (headless Octave run)

Headless adaptation of `co2lab/co2lab-spillpoint/examples/firstTrappingExample.m`,
run on 2026-08-19 with GNU Octave 8.4.

- `run_first_trapping.m` — the driver. It uses the edge-based method in
  `trapAnalysis`, because the cell-based method needs `matlab_bgl` MEX
  binaries, which Octave cannot load. It uses 2-D maps, because the gnuplot
  toolkit cannot draw 3-D quad patches.
- `run.log` — the run log.
- Figures: depth map with trap cells, traps with spill paths, trap volumes.

Result: 13 structural traps; total trapping capacity 5.099e+06 m3, which is
40.79 % of the total pore volume (1.25e+07 m3).
