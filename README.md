Experiments attempting extending the Winslow high-order mesh smoother of
Fortunato & Persson (*JCP* 2016) with a `W`-modification intended to
decouple the smoothed physical mesh from the shape of the reference
(computational) mesh.

## `w1.m`, `w2.m`, `w2_3d.m`, `w2_skewed.m` — the original Fortunato–Persson paradigm

- `init.m` — sets the MATLAB path; run before any demo.
- `elliptic_smoothing.m` — the 2016 Fortunato–Persson Winslow smoother

These scripts reproduce the unmodified 2016 formulation: given a
valid (untangled) **reference mesh**, prescribe a boundary perturbation
that would produce a heavily tangled initial **physical mesh**, then
run the Winslow iteration until all elements are untangled and
well-shaped. All four live in `demos/original-formulation/` and save
figures to `demos/original-formulation/figures/<script>/`.

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
(`mshsquare(n+1, n+1)`) at `porder = 4`. The boundary on side 4 (`x = 0`)
is displaced by `0.5 * sin(π y)`, cutting roughly five element layers
into the domain and producing heavy interior tangling in the initial
configuration. Three figures are saved to `demos/original-formulation/figures/w2/`:
`mesh_reference.png` (the clean reference grid), `mesh_tangled.png`
(the initial physical configuration, with many inverted elements near
side~4), and `mesh_smoothed.png` (the Winslow result, with all
elements valid). This is the canonical
reference-mesh-is-good / physical-mesh-starts-tangled test and is the
starting point for all the subsequent experiments in this repo.

### `w2_3d.m` — 3D analogue of `w2.m`

Same paradigm in three dimensions. Reference is a structured
tetrahedral mesh of the unit cube (`mshcube(n+1, n+1, n+1)` with
`n = 6`, giving a `7 × 7 × 7` vertex grid and 1296 tets) at
`porder = 4`. The boundary on side~1 (`x = 0`) is displaced inward
by `0.4 * sin(π y) * sin(π z)`, producing a severely tangled initial
configuration near that face (`min I ≈ -7`, where `I` is the scaled
Jacobian of Fortunato–Persson). After the Winslow iteration all
elements are valid (`min I ≈ +0.7`). The script saves three mesh
views to `demos/original-formulation/figures/w2_3d/` (reference, tangled, smoothed), all
rendered from a viewpoint that shows the perturbed `x = 0` face,
along with three per-element scaled-Jacobian histograms (one per
configuration) demonstrating that the negative tail of `I` in the
tangled mesh is eliminated after smoothing.

### `w2_skewed.m` — `w2.m` with a slightly skewed reference

A minimal variant of `w2.m` that serves as the main test case for the
upcoming `W`-modification experiments. The setup is identical to
`w2.m`—same `11 × 11` grid, same `porder = 4`, same side-4 boundary
perturbation—with one exception: one interior linear vertex of the
reference is manually displaced toward a neighbour, creating a small
cluster of sliver-shaped elements in the upper-right of the reference.
Unmodified Fortunato–Persson still converges to a valid physical
mesh, but the sliver cluster from the reference persists as the
worst-shaped element in the smoothed output (`max η ≈ 7.2` on the
sliver, vs `1.0` on the rest of the mesh). The script saves the same
three mesh figures as `w2.m` plus two histogram triptychs in the same
scaled-Jacobian `I` and Liu–Joe shape-distortion `η` metrics.
Because the scaled Jacobian `I` is blind to straight-sided distortion
(`I = 1` identically on the reference regardless of element shape),
`η` is the metric that actually sees the sliver. This inheritance of
reference distortion by the smoothed output is exactly the failure
mode that motivates the `W`-modification.

## ADD NEXT BLOCK OF EXPERIMENTS.
