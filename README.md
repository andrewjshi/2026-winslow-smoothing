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

## `w2_W.m`, `w2_skewed_W.m` — the `W`-modified formulation

The next batch of experiments modifies the split-form metric
`g^ij → g̃^ij = W_h g^ij W_h^T`, where `W_h` is a continuous nodal
field built from per-element target Jacobians mapping an ideal
equilateral triangle to each reference element. The smoother's
Picard structure is left unchanged; only the right-hand side
assembly sees the modification.

**Auxiliary functions** (live at the project root):

- `elliptic_smoothing_W.m` — a copy of `elliptic_smoothing.m` with
  an `isfield(msh, 'W_h')` guard inside `sm1alpha2d`/`sm12d`. When
  `msh.W_h` is not set, the code runs the 2016 formulation
  bit-for-bit; otherwise it substitutes `g̃^ij` as described above.
- `compute_WK_simplex.m` — builds the per-element constant `W_K`
  from the straight-sided reference corners and an ideal element
  (default: equilateral simplex).
- `compute_Wh_nodal.m` — projects `W_K` onto a nodal `W_h` field in
  `V_h^p` (default) or `V_h^1`, using either a full consistent mass
  matrix (default, matching the original `α`-projection) or a
  lumped-mass (`Eq. 9`-style area-weighted) simplification.

**Drivers** (live in `demos/w-formulation/` and save figures to
`demos/w-formulation/figures/<script>/`):

- `w2_W.m` — `w2.m` run through the `W`-modified solver. Starting
  point closest to the 2016 formulation: `V_h^p + full consistent
  mass`. First result is that this produces two interior elements
  with `min J < 0` on a problem where the unmodified solver
  produces a valid mesh.
- `w2_skewed_W.m` — `w2_skewed.m` run through the same `W`-modified
  solver. The main test: does the `W`-modification cancel the
  sliver distortion the unmodified solver inherits from the
  reference?
- `w_regression.m` — bit-for-bit regression. Verifies that
  `W_h = I` produces the exact same output as `elliptic_smoothing`.
  Run after any change to the W machinery.
