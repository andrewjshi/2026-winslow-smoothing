# 2026-winslow-smoothing

Experiments extending the Winslow-Picard high-order mesh smoother of
Fortunato & Persson (*JCP* 2016) with a `W`-modification intended to
decouple the smoothed physical mesh from the shape of the reference
(computational) mesh.

## Repo layout

- `dgmatlab/` — Per-Olof Persson's DG toolkit (shipped unchanged, minus a
  pure-MATLAB replacement for `dgmass` since the bundled MEX binaries do
  not run on Apple Silicon; the replacement is `dgmass_ml.m` at the repo
  root).
- `elliptic_smoothing.m` — the 2016 Fortunato–Persson Winslow-Picard
  smoother (minor patch at the project root to call `dgmass_ml` instead
  of the MEX `dgmass`).
- `elliptic_smoothing_W.m` — copy of the above with the `W`-modification
  guarded behind `isfield(msh, 'W_h')`; falls back to the unmodified
  code when `msh.W_h` is absent.
- `compute_WK_simplex.m` — per-element `W_K` from straight-sided corners
  to an ideal (equilateral) simplex.
- `compute_Wh_nodal.m` — lumped L² projection of `W_K` onto a continuous
  nodal field `W_h` (P¹ or full V_h^p, with area / η / η²-style
  weightings).
- `demos/` — MATLAB scripts documenting the experiments below.
- `notes/winslow-notes.tex` — the companion writeup.
- `init.m` — sets the MATLAB path; run before any demo.

## Running

```matlab
run('init.m')
cd demos
w2
```

## `w1.m` and `w2.m` — the original Fortunato–Persson paradigm

These two scripts reproduce the unmodified 2016 formulation: given a
valid (untangled) **reference mesh**, prescribe a boundary perturbation
that would produce a heavily tangled initial **physical mesh**, then
run the Winslow–Picard iteration until all elements are untangled and
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
side~4), and `mesh_smoothed.png` (the Winslow-iterated result, with
all elements valid). This is the canonical
reference-mesh-is-good / physical-mesh-starts-tangled test and is the
starting point for all the subsequent experiments in this repo.
