% stadium_sliding.m - Boundary-sliding configurations on the convex stadium.
%
% Pipeline: same stadium mesh as stadium_preview.m, then run four
% configurations of the shape-opt + Winslow handoff and compare quality:
%
%   A. Pinned shape-opt, pinned Winslow                   (baseline)
%   B. Sliding shape-opt, pinned Winslow                  (with pullback)
%   C. Pinned shape-opt, sliding Winslow                  (curve sliding)
%   D. Sliding shape-opt, sliding Winslow                 (both)
%
% The stadium boundary is C^1 (no real corners), so we pin a single
% reference vertex bnd_v(1) on both sides to lock the parametrization.

rng(3);
porder = 4;
a = 1.0;
b = 0.5;
n_bnd = 80;
n_int = 100;

% ---------- build stadium mesh (matches stadium_preview.m) ---------------
perim_straight = 2*a;
perim_semi     = pi*b;
P_stadium      = 2*perim_straight + 2*perim_semi;
s_vals = P_stadium * (0:n_bnd-1)' / n_bnd;
bnd = zeros(n_bnd, 2);
for i = 1:n_bnd
    s = s_vals(i);
    if s < perim_straight
        bnd(i,:) = [-a + s, -b];
    elseif s < perim_straight + perim_semi
        t = (s - perim_straight)/b - pi/2;
        bnd(i,:) = [a + b*cos(t), b*sin(t)];
    elseif s < 2*perim_straight + perim_semi
        bnd(i,:) = [a - (s - perim_straight - perim_semi), b];
    else
        t = (s - 2*perim_straight - perim_semi)/b + pi/2;
        bnd(i,:) = [-a + b*cos(t), b*sin(t)];
    end
end

stadium_fd = @(p) sqrt(max(abs(p(:,1)) - a, 0).^2 + p(:,2).^2) - b;

x_min = -(a + b) - 0.02;  x_max = (a + b) + 0.02;
y_min = -b - 0.02;         y_max = b + 0.02;
int_pts = zeros(0, 2);
while size(int_pts, 1) < n_int
    pt = [x_min + (x_max - x_min) * rand(), ...
          y_min + (y_max - y_min) * rand()];
    if stadium_fd(pt) < -0.02
        int_pts(end+1, :) = pt; %#ok<AGROW>
    end
end

all_pts = [bnd; int_pts];
n_b = size(bnd, 1);
ce = [(1:n_b)', [(2:n_b)'; 1]];
dt = delaunayTriangulation(all_pts, ce);
in = isInterior(dt);
tri = dt.ConnectivityList(in, :);

bndexpr = {'true'};
msh = ml2msh(all_pts, tri, bndexpr, stadium_fd, {});
msh = mshcurved(msh, 1);
msh = nodealloc(msh, porder);

figdir = 'figures/stadium_sliding';
if ~exist(figdir, 'dir'), mkdir(figdir); end

% ---------- Tutte embedding onto unit disk -------------------------------
msh_t = tutte_embedding(msh, 'TargetShape', 'circle');
msh_t = nodealloc(msh_t, porder);

% Identify the boundary loop as ordered by tutte_embedding so we can
% choose pin points spread around the perimeter. The first vertex of
% the CCW-oriented loop is at theta=0 on the disk; subsequent vertices
% follow CCW.
bnd_v = extract_bnd_loop(msh_t);
if signed_area_local(msh_t.p(:, bnd_v)) < 0
    bnd_v = flip(bnd_v);
end
th_loop = atan2(msh_t.p(2, bnd_v), msh_t.p(1, bnd_v));
[~, k0] = min(abs(mod(th_loop, 2*pi)));
bnd_v = bnd_v([k0:end, 1:k0-1]);

% Sliding modes need pinned anchors to keep the elliptic PDE
% well-conditioned and the parametrization stable. With a curved C^1
% boundary, per-CG-node projection on too-loose pinning scrambles the
% high-order DG nodes; pin 8 evenly-spaced linear vertices to anchor
% the boundary.
nbnd_loop = numel(bnd_v);
n_pins = 8;
pin_idx = round((0:n_pins-1) * nbnd_loop / n_pins) + 1;
pin_idx = unique(min(pin_idx, nbnd_loop));
pin_v   = bnd_v(pin_idx);
pin_xy  = msh.p(:, pin_v)';
ref_v   = pin_v(1);
fprintf('  pinned anchor vertices: idx %s\n', mat2str(pin_v(:)'));
for kk = 1:numel(pin_v)
    fprintf('    (%g, %g)\n', pin_xy(kk,1), pin_xy(kk,2));
end

stadium_proj = @(p) stadium_project(p, a, b);

plot_mesh(msh,    fullfile(figdir, 'mesh_clean.png'),         [900 400]);
plot_mesh(msh_t,  fullfile(figdir, 'mesh_tutte_circle.png'),  [600 600]);

% ---------- A. Pinned shape-opt + pinned Winslow -------------------------
fprintf('\n=== A: pinned / pinned ===\n');
[msh_rA, infoA] = optimize_shape(msh_t);
msh_rA = nodealloc(msh_rA, porder);
fprintf('  shape-opt: F0=%.4e -> F=%.4e, iters=%d\n', ...
    infoA.F_initial, infoA.F_final, infoA.iterations);
plot_mesh(msh_rA, fullfile(figdir, 'A_optimized.png'), [600 600]);

is_bnd_dg_A = boundary_dg_mask(msh_rA);
curvep1 = msh_rA.p1;
curvep1(is_bnd_dg_A) = msh.p1(is_bnd_dg_A);
msh_wA = elliptic_smoothing(msh_rA, curvep1, false);
plot_mesh(msh_wA, fullfile(figdir, 'A_winslow.png'),  [900 400]);

% ---------- B. Sliding shape-opt + pinned Winslow ------------------------
fprintf('\n=== B: slide / pinned ===\n');
[msh_rB, infoB] = optimize_shape(msh_t, ...
    'SlideBoundary', 'circle', 'PinBoundaryNodes', pin_v);
msh_rB = nodealloc(msh_rB, porder);
fprintf('  shape-opt: F0=%.4e -> F=%.4e, iters=%d\n', ...
    infoB.F_initial, infoB.F_final, infoB.iterations);
plot_mesh(msh_rB, fullfile(figdir, 'B_optimized.png'), [600 600]);

% Derive new physical boundary positions via proportional-arc-length pullback.
msh_phys_B = derive_physical_from_slid_circle( ...
    msh_rB, msh, porder, ref_v, stadium_fd);
plot_mesh(msh_phys_B, fullfile(figdir, 'B_physical_target.png'), [900 400]);

is_bnd_dg_B = boundary_dg_mask(msh_rB);
curvep1 = msh_rB.p1;
curvep1(is_bnd_dg_B) = msh_phys_B.p1(is_bnd_dg_B);
msh_wB = elliptic_smoothing(msh_rB, curvep1, false);
plot_mesh(msh_wB, fullfile(figdir, 'B_winslow.png'), [900 400]);

% ---------- C. Pinned shape-opt + sliding Winslow ------------------------
fprintf('\n=== C: pinned / slide ===\n');
spec_C.project_fn = stadium_proj;
spec_C.pinned_xy  = pin_xy;
curvep1 = msh_rA.p1;
curvep1(is_bnd_dg_A) = msh.p1(is_bnd_dg_A);
msh_wC = elliptic_smoothing_slide(msh_rA, curvep1, false, spec_C);
plot_mesh(msh_wC, fullfile(figdir, 'C_winslow.png'), [900 400]);

% ---------- D. Sliding shape-opt + sliding Winslow -----------------------
fprintf('\n=== D: slide / slide ===\n');
spec_D.project_fn = stadium_proj;
spec_D.pinned_xy  = pin_xy;
curvep1 = msh_rB.p1;
curvep1(is_bnd_dg_B) = msh_phys_B.p1(is_bnd_dg_B);
msh_wD = elliptic_smoothing_slide(msh_rB, curvep1, false, spec_D);
plot_mesh(msh_wD, fullfile(figdir, 'D_winslow.png'), [900 400]);

% ---------- Quality summary ---------------------------------------------
fprintf('\n=== quality summary (Winslow output, per-element eta) ===\n');
fprintf('  %-3s  %10s  %10s  %12s  %10s\n', ...
    'cfg', 'mean_eta', 'max_eta', 'sum_eta^2', 'min_I');
report_quality('A', msh_wA);
report_quality('B', msh_wB);
report_quality('C', msh_wC);
report_quality('D', msh_wD);

fprintf('\nSaved figures to %s/\n', figdir);


% =========================================================================
% Helpers
% =========================================================================

function q = stadium_project(p, a, b)
% Closest point on the stadium boundary to p (1x2 row vector).
if p(1) > a
    d = p - [a, 0];
    q = [a, 0] + b * d / norm(d);
elseif p(1) < -a
    d = p - [-a, 0];
    q = [-a, 0] + b * d / norm(d);
else
    if p(2) >= 0
        q = [p(1),  b];
    else
        q = [p(1), -b];
    end
end
end


function msh_phys_new = derive_physical_from_slid_circle(msh_r, msh_phys, porder, ref_v, fd)
% Given a shape-opt output `msh_r` with boundary linear vertices on the
% unit circle (with ref_v pinned at theta=0), derive new physical boundary
% positions on `msh_phys`'s domain via proportional-arc-length pullback.
% Re-curves the result via mshcurved using `fd` so high-order DG nodes
% land on the analytic boundary curve.
bnd_v = extract_bnd_loop(msh_r);
if signed_area_local(msh_r.p(:, bnd_v)) < 0
    bnd_v = flip(bnd_v);
end
i_ref = find(bnd_v == ref_v, 1);
if isempty(i_ref)
    error('derive_physical_from_slid_circle: ref_v not found on boundary.');
end
bnd_v = bnd_v([i_ref:end, 1:i_ref-1]);
nbnd = numel(bnd_v);

edge_lens = vecnorm( ...
    msh_phys.p(:, bnd_v([2:end, 1])) - msh_phys.p(:, bnd_v), 2, 1)';
cum_arc = [0; cumsum(edge_lens)];
P = cum_arc(end);

slid_pos = msh_r.p(:, bnd_v);
th = atan2(slid_pos(2,:), slid_pos(1,:))';
th_rel = mod(th - th(1), 2*pi);
th_rel(1) = 0;

s_new = th_rel * P / (2*pi);

new_bnd_pos = zeros(2, nbnd);
for i = 1:nbnd
    s = s_new(i);
    k = find(cum_arc(1:end-1) <= s & s < cum_arc(2:end), 1);
    if isempty(k), k = nbnd; end
    denom = cum_arc(k+1) - cum_arc(k);
    if denom < 1e-12
        frac = 0;
    else
        frac = (s - cum_arc(k)) / denom;
    end
    p1 = msh_phys.p(:, bnd_v(k));
    p2 = msh_phys.p(:, bnd_v(mod(k, nbnd) + 1));
    new_bnd_pos(:, i) = p1 + frac * (p2 - p1);
end

msh_phys_new = msh_phys;
msh_phys_new.p(:, bnd_v) = new_bnd_pos;
if nargin >= 5 && ~isempty(fd)
    msh_phys_new.fd = fd;
    msh_phys_new.fdargs = {};
    msh_phys_new = mshcurved(msh_phys_new, 1);
end
msh_phys_new = nodealloc(msh_phys_new, porder);
end


function bnd_v = extract_bnd_loop(msh)
t = double(msh.t) + 1;
t2t = double(msh.t2t);
edges = [];
for it = 1:size(t, 2)
    for j = 1:size(t2t, 1)
        if t2t(j, it) < 0
            v = t(:, it); v(j) = [];
            edges = [edges; v(:)']; %#ok<AGROW>
        end
    end
end
nedges = size(edges, 1);
bnd_v = zeros(nedges, 1);
bnd_v(1) = edges(1, 1);
curr = edges(1, 2);
used = false(nedges, 1);
used(1) = true;
for k = 2:nedges
    bnd_v(k) = curr;
    for ie = 1:nedges
        if used(ie), continue; end
        if edges(ie, 1) == curr
            curr = edges(ie, 2); used(ie) = true; break;
        elseif edges(ie, 2) == curr
            curr = edges(ie, 1); used(ie) = true; break;
        end
    end
end
end


function s = signed_area_local(p)
n = size(p, 2);
s = 0;
for k = 1:n
    kn = mod(k, n) + 1;
    s = s + p(1, k) * p(2, kn) - p(1, kn) * p(2, k);
end
s = s / 2;
end


function eta = tri_eta_linear(p, t)
nt = size(t, 2);
eta = zeros(nt, 1);
A1 = 4/3; A2 = -2/3; A3 = 4/3;
det_Minv = 2/sqrt(3);
for it = 1:nt
    p1 = p(:, t(1, it));  p2 = p(:, t(2, it));  p3 = p(:, t(3, it));
    a_ = p2(1) - p1(1);  b_ = p3(1) - p1(1);
    c_ = p2(2) - p1(2);  d_ = p3(2) - p1(2);
    detE = a_ * d_ - b_ * c_;
    if detE <= 0
        eta(it) = Inf; continue;
    end
    Q = (a_^2 + c_^2) * A1 + 2 * (a_*b_ + c_*d_) * A2 + (b_^2 + d_^2) * A3;
    D = 2 * det_Minv * detE;
    eta(it) = Q / D;
end
end


function I = tri_I_linear(p, t)
% Per-element invertibility marker: detE>0 -> I=1, else I=0.
nt = size(t, 2);
I = zeros(nt, 1);
for it = 1:nt
    p1 = p(:, t(1, it));  p2 = p(:, t(2, it));  p3 = p(:, t(3, it));
    a_ = p2(1) - p1(1);  b_ = p3(1) - p1(1);
    c_ = p2(2) - p1(2);  d_ = p3(2) - p1(2);
    if a_*d_ - b_*c_ > 0, I(it) = 1; end
end
end


function report_quality(label, msh)
eta = tri_eta_linear(msh.p, double(msh.t) + 1);
I   = tri_I_linear(msh.p, double(msh.t) + 1);
fprintf('  %-3s  %10.4f  %10.4f  %12.2f  %10.0f\n', ...
    label, mean(eta), max(eta), sum(eta.^2), min(I));
end


function plot_mesh(msh, filename, figsize)
figure; clf;
dgmeshplot_curved(msh, 4, 1, 0);
set(findobj(gca, 'Type', 'line'), 'MarkerSize', 1);
axis equal; axis off;
set(gcf, 'Position', [100 100 figsize(1) figsize(2)]);
exportgraphics(gcf, filename, 'Resolution', 200);
close;
end


function m = boundary_dg_mask(msh)
data = dginit(msh);
[~, edg] = getbndnodes(msh, data);
[ns, d, nt] = size(msh.p1);
m2 = false(ns, nt);
m2(edg) = true;
m = repmat(reshape(m2, ns, 1, nt), 1, d, 1);
end
