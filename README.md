Experiments attempting extending the Winslow high-order mesh smoother of
Fortunato & Persson (*JCP* 2016) with a `W`-modification intended to
decouple the smoothed physical mesh from the shape of the reference
(computational) mesh.

## `w1.m` and `w2.m` — the original Fortunato–Persson paradigm

- `init.m` — sets the MATLAB path; run before any demo.
- `elliptic_smoothing.m` — the 2016 Fortunato–Persson Winslow smoother
  (minor patch at the project root to call `dgmass_ml` instead of the
  MEX `dgmass`, since the shipped MEX binaries do not run on Apple
  Silicon).

These two scripts reproduce the unmodified 2016 formulation: given a
valid (untangled) **reference mesh**, prescribe a boundary perturbation
that would produce a heavily tangled initial **physical mesh**, then
run the Winslow iteration until all elements are untangled and
well-shaped.

### `w1.m` — curved boundary, unit disk

Reference mesh from `mshcircle(1)` with `porder = 4`. The linear
boundary is promoted to the true circular boundary via `mshcurved`, so
the initial physical high-order mesh has curved edges that fold onto
themselves near the boundary. `elliptic_smoothing` untangles them and
returns a valid curved mesh on the disk. This is the simplest
demonstration: the reference is a clean triangulation of the disk;
the physical target is the same disk with curved high-order edges
prescribed on the boundary.

### `w2.m` — structured square, violent boundary indent

Reference is a structured `11 × 11` uniform mesh of the unit square
(`mshsquare(n+1, n+1)`) at `porder = 4`. The boundary on side~4 (`x = 0`)
is displaced by `0.5 * sin(π y)`, cutting roughly five element layers
into the domain and producing heavy interior tangling in the initial
configuration. Three figures are saved to `demos/figures/w2/`:
`mesh_reference.png` (the clean reference grid), `mesh_tangled.png`
(the initial physical configuration, with many inverted elements near
side~4), and `mesh_smoothed.png` (the Winslow result, with all
elements valid). This is the canonical
reference-mesh-is-good / physical-mesh-starts-tangled test and is the
starting point for all the subsequent experiments in this repo.

## ADD NEXT BLOCK OF EXPERIMENTS.
