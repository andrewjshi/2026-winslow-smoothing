% l_shape_tangled.m - Tangled L-shape input through the Tutte pipeline,
% tested against both unit-square and unit-circle computational domains.
%
% Pipeline:
%   1. Build unstructured Delaunay mesh of an L-shape.
%   2. Tangle interior nodes violently (both linear and high-order) so
%      they can end up outside the physical domain.
%   3. For each target ('square', 'circle'):
%        a. Tutte embed on the target.
%        b. Shape-optimise.
%        c. nodealloc to high order.
%      [FP hand-off to elliptic_smoothing is deferred for now --- we just
%       want to see the reference-mesh construction here.]

rng(17);
porder = 4;

% ---------------------------------------------------------------------
% Step 1: build a clean unstructured L-shape mesh.
% L = unit square minus [0.5, 1] x [0.5, 1].
% Corners (CCW): (0,0), (1,0), (1,0.5), (0.5,0.5), (0.5,1), (0,1).
% ---------------------------------------------------------------------
n_side     = 5;       % boundary segments per L-side
n_interior = 35;

% Boundary points on each of the 6 sides (open intervals, no corners yet)
segs = [
    0   0   1   0;    % bottom (x from 0 to 1, y=0)
    1   0   1   0.5;
    1   0.5 0.5 0.5;
    0.5 0.5 0.5 1;
    0.5 1   0   1;
    0   1   0   0];
bnd_pts = [];
for s = 1:size(segs, 1)
    p1 = segs(s, 1:2); p2 = segs(s, 3:4);
    for k = 0:n_side-1
        t = k / n_side;
        bnd_pts = [bnd_pts; p1 + t * (p2 - p1)]; %#ok<AGROW>
    end
end

% Interior points: random in the L (reject the cutout)
int_pts = zeros(0, 2);
while size(int_pts, 1) < n_interior
    pt = 0.05 + 0.9 * rand(1, 2);
    if pt(1) > 0.5 && pt(2) > 0.5
        continue;   % in the cut-out
    end
    int_pts = [int_pts; pt]; %#ok<AGROW>
end

% Anchor points near each L corner to guarantee 3-connectedness (each
% boundary corner must have at least one interior neighbour in the
% Delaunay).
corner_anchors = [0.08 0.08; 0.92 0.08; 0.92 0.42; 0.42 0.42; 0.42 0.92; 0.08 0.92];
int_pts = [int_pts; corner_anchors];

all_pts = [bnd_pts; int_pts];

% Constrained Delaunay respecting the L's boundary loop.
n_bnd = size(bnd_pts, 1);
ce = [(1:n_bnd)', [(2:n_bnd)'; 1]];
dt = delaunayTriangulation(all_pts, ce);
in = isInterior(dt);
tri = dt.ConnectivityList(in, :);

bndexpr = {
    'all(p(:,2) < 1e-6)'                                         % bottom
    'all(p(:,1) > 1 - 1e-6)'                                     % lower right
    'all(p(:,2) > 0.5 - 1e-6 & p(:,2) < 0.5 + 1e-6)'             % inner top
    'all(p(:,1) > 0.5 - 1e-6 & p(:,1) < 0.5 + 1e-6)'             % inner right
    'all(p(:,2) > 1 - 1e-6)'                                     % upper top
    'all(p(:,1) < 1e-6)'                                         % left
};
msh_clean = ml2msh(all_pts, tri, bndexpr);
msh_clean = mshcurved(msh_clean, []);
msh_clean = nodealloc(msh_clean, porder);

figdir = 'figures/l_shape_tangled';
if ~exist(figdir, 'dir'), mkdir(figdir); end

% --- 6 L-shape domain corners (as linear-vertex indices); used to
% track the corners through Tutte and shape-opt in figure panels.
lshape_corners = [0 0; 1 0; 1 0.5; 0.5 0.5; 0.5 1; 0 1];
n_lc = size(lshape_corners, 1);
lshape_corner_idx = zeros(n_lc, 1);
for k = 1:n_lc
    [~, lshape_corner_idx(k)] = min(vecnorm( ...
        msh_clean.p' - lshape_corners(k,:), 2, 2));
end

fprintf('clean L-shape mesh:\n'); print_quality(msh_clean);
plot_mesh_verts(msh_clean, lshape_corner_idx, ...
    fullfile(figdir, 'mesh_clean.png'));

% --- 3-connectedness diagnostic: find any pair whose removal disconnects
[ok, pair] = is_3_connected_diagnostic(msh_clean);
if ok
    fprintf('graph is 3-connected.\n');
else
    fprintf('graph NOT 3-connected; removing vertices %d, %d disconnects it.\n', ...
            pair(1), pair(2));
    fprintf('  positions: (%.3f, %.3f) and (%.3f, %.3f)\n', ...
            msh_clean.p(1, pair(1)), msh_clean.p(2, pair(1)), ...
            msh_clean.p(1, pair(2)), msh_clean.p(2, pair(2)));
end

% ---------------------------------------------------------------------
% Step 2: tangle the mesh violently.
% Move INTERIOR linear vertices by large random displacements (may land
% outside the L). Also displace high-order nodes. Boundary preserved.
% ---------------------------------------------------------------------
msh_tangled = msh_clean;

is_bnd_lin = boundary_linear_mask(msh_tangled);
amp = 0.25;
d_lin = amp * (rand(size(msh_tangled.p)) - 0.5) * 2;   % ~ Uniform[-amp, amp]
d_lin(:, is_bnd_lin) = 0;
msh_tangled.p = msh_tangled.p + d_lin;
msh_tangled = nodealloc(msh_tangled, porder);

% Additional high-order jitter on interior DG nodes only
is_bnd_dg = boundary_dg_mask(msh_tangled);
d_ho = 0.08 * (rand(size(msh_tangled.p1)) - 0.5) * 2;
d_ho(is_bnd_dg) = 0;
msh_tangled.p1 = msh_tangled.p1 + d_ho;

fprintf('\ntangled input:\n'); print_quality(msh_tangled);
plot_mesh(msh_tangled, fullfile(figdir, 'mesh_tangled.png'));

% ---------------------------------------------------------------------
% Step 3: for each (target in {square, circle}) x (slide in {pinned,
% slide}), run Tutte + shape-opt + Winslow hand-off, and record quality
% metrics at both the reference-mesh stage (post shape-opt) and the
% physical-mesh stage (post Winslow).
% ---------------------------------------------------------------------
summary = struct('target', {}, 'mode', {}, 'stage', {}, ...
    'mean_eta', {}, 'max_eta', {}, 'mean_eta2', {}, 'max_eta2', {});

square_corners = [0 1 1 0; 0 0 1 1];

for target_shape = {'square', 'circle'}
    ts = target_shape{1};
    fprintf('\n=== target: %s ===\n', ts);

    msh_t = tutte_embedding(msh_tangled, 'TargetShape', ts);
    msh_t = nodealloc(msh_t, porder);
    fprintf('after Tutte:\n'); print_quality(msh_t);
    plot_mesh_verts(msh_t, lshape_corner_idx, ...
        fullfile(figdir, sprintf('mesh_tutte_%s.png', ts)));

    for mode_cell = {'pinned', 'slide'}
        mode = mode_cell{1};
        fprintf('\n--- %s / %s ---\n', ts, mode);

        if strcmp(mode, 'pinned')
            [msh_r, info] = optimize_shape(msh_t);
        elseif strcmp(ts, 'square')
            [msh_r, info] = optimize_shape(msh_t, ...
                'SlideBoundary', 'polygon', 'TargetCorners', square_corners);
        else
            % Circle slide with L-shape corners pinned on the circle at
            % their proportional-arc positions.
            [msh_r, info] = optimize_shape(msh_t, ...
                'SlideBoundary', 'circle', ...
                'PinBoundaryNodes', lshape_corner_idx);
        end
        msh_r = nodealloc(msh_r, porder);

        fprintf('after shape-opt (%s):\n', mode); print_quality(msh_r);
        fprintf('  F_init = %.4e, F_final = %.4e, iters = %d\n', ...
                info.F_initial, info.F_final, info.iterations);
        if strcmp(mode, 'pinned')
            plot_mesh_verts(msh_r, lshape_corner_idx, fullfile(figdir, ...
                sprintf('mesh_optimized_%s.png', ts)));
        else
            plot_mesh_verts(msh_r, lshape_corner_idx, fullfile(figdir, ...
                sprintf('mesh_optimized_%s_slide.png', ts)));
        end
        summary(end+1) = pack_row(ts, mode, 'shape-opt', msh_r);

        % Winslow hand-off.
        %
        % Pinned mode: Dirichlet data is the tangled input's L-shape
        % boundary, i.e., the uniform arc-length sampling of the original L.
        %
        % Slide mode (circle only): derive NEW L-shape boundary positions
        % that correspond to the slid circle positions via proportional
        % arc-length mapping. This keeps the reference-physical
        % correspondence consistent so pinned Winslow converges.
        is_bnd_dg_ref = boundary_dg_mask(msh_r);
        curvep1 = msh_r.p1;
        if strcmp(mode, 'slide') && strcmp(ts, 'circle')
            msh_phys_new = derive_physical_from_slid_circle( ...
                msh_r, msh_tangled, porder, lshape_corner_idx(1));
            curvep1(is_bnd_dg_ref) = msh_phys_new.p1(is_bnd_dg_ref);
        else
            curvep1(is_bnd_dg_ref) = msh_tangled.p1(is_bnd_dg_ref);
        end
        try
            msh_w = elliptic_smoothing(msh_r, curvep1, false);
            fprintf('after Winslow (%s):\n', mode); print_quality(msh_w);
            [~, I_w] = quality(msh_w);
            inv_tris = find(I_w <= 0);
            if ~isempty(inv_tris)
                fprintf('  inverted element indices: %s\n', mat2str(inv_tris(:)'));
            end
            if strcmp(mode, 'pinned')
                plot_mesh_hl(msh_w, inv_tris, fullfile(figdir, ...
                    sprintf('mesh_winslow_%s.png', ts)));
            else
                plot_mesh_hl(msh_w, inv_tris, fullfile(figdir, ...
                    sprintf('mesh_winslow_%s_slide.png', ts)));
            end
            summary(end+1) = pack_row(ts, mode, 'winslow', msh_w);
        catch err
            fprintf('Winslow FAILED (%s / %s): %s\n', ts, mode, err.message);
            summary(end+1) = struct('target', ts, 'mode', mode, ...
                'stage', 'winslow-FAIL', ...
                'mean_eta', NaN, 'max_eta', NaN, ...
                'mean_eta2', NaN, 'max_eta2', NaN);
        end
    end
end

% --- Summary table ---
fprintf('\n=== quality summary (per element) ===\n');
fprintf('  %-8s  %-6s  %-10s  %8s  %8s  %10s  %10s\n', ...
    'target', 'mode', 'stage', 'mean_eta', 'max_eta', 'mean_eta2', 'max_eta2');
for s = summary
    fprintf('  %-8s  %-6s  %-10s  %8.4f  %8.4f  %10.4f  %10.4f\n', ...
        s.target, s.mode, s.stage, s.mean_eta, s.max_eta, s.mean_eta2, s.max_eta2);
end

fprintf('\nSaved figures to %s/\n', figdir);


% ====================================================================
function print_quality(msh)
[eta, I] = quality(msh);
fprintf('  nelem = %d\n', numel(eta));
fprintf('  eta   : min = %.4g, mean = %.4g, max = %.4g\n', ...
        min(eta), mean(eta), max(eta));
fprintf('  I     : min = %.4g  (I<=0 => inverted)\n', min(I));
end


function msh_phys_new = derive_physical_from_slid_circle(msh_r, msh_tangled, porder, ref_v)
% Given a shape-opt output `msh_r` with boundary linear vertices on the
% unit circle (some slid), derive new L-shape boundary positions via
% proportional arc-length mapping and return a physical mesh with those
% new boundary positions. `ref_v` is a linear-vertex index of a PINNED
% boundary vertex (e.g., an L-shape corner) used as the arc-length
% reference. Using a pinned vertex is essential: its angle on the circle
% is unchanged by sliding, so the cum_arc reference frame doesn't drift.

bnd_v = extract_bnd_loop(msh_r);
if signed_area_local(msh_r.p(:, bnd_v)) < 0
    bnd_v = flip(bnd_v);
end
% Rotate bnd_v so that ref_v is first
i_ref = find(bnd_v == ref_v, 1);
if isempty(i_ref)
    error('derive_physical_from_slid_circle: ref_v not found on boundary.');
end
bnd_v = bnd_v([i_ref:end, 1:i_ref-1]);
nbnd = numel(bnd_v);

% Cumulative arc length on the ORIGINAL L-shape perimeter (waypoint at
% each boundary linear vertex, interpolated linearly between).
edge_lens = vecnorm( ...
    msh_tangled.p(:, bnd_v([2:end, 1])) - msh_tangled.p(:, bnd_v), 2, 1)';
cum_arc = [0; cumsum(edge_lens)];
P = cum_arc(end);

% Slid angles on circle, relative to bnd_v(1) so we don't pick up an
% arbitrary offset. The angle of bnd_v(1) is already 0 by construction of
% the proportional-arc Tutte embedding, but we subtract anyway.
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
    p1 = msh_tangled.p(:, bnd_v(k));
    p2 = msh_tangled.p(:, bnd_v(mod(k, nbnd) + 1));
    new_bnd_pos(:, i) = p1 + frac * (p2 - p1);
end

msh_phys_new = msh_tangled;
msh_phys_new.p(:, bnd_v) = new_bnd_pos;
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
            edges = [edges; v(:)'];  %#ok<AGROW>
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


function row = pack_row(target, mode, stage, msh)
[eta, ~] = quality(msh);
row = struct('target', target, 'mode', mode, 'stage', stage, ...
    'mean_eta', mean(eta), 'max_eta', max(eta), ...
    'mean_eta2', mean(eta .^ 2), 'max_eta2', max(eta .^ 2));
end


function val = mean_eta2(msh)
[eta, ~] = quality(msh);
val = mean(eta .^ 2);
end


function [eta, I] = quality(msh)
d = size(msh.p, 1);
data = dginit(msh);
p1x = permute(msh.p1(:,1,:), [1,3,2]);
p1y = permute(msh.p1(:,2,:), [1,3,2]);
phiX = permute(data.gfs(:,2,:), [1,3,2]);
phiY = permute(data.gfs(:,3,:), [1,3,2]);
xX = phiX'*p1x; xY = phiY'*p1x;
yX = phiX'*p1y; yY = phiY'*p1y;
detJ = xX.*yY - xY.*yX;
Jfro2 = xX.^2 + xY.^2 + yX.^2 + yY.^2;
etap = Jfro2 ./ (d * max(abs(detJ), eps).^(2/d));
eta = max(etap, [], 1)';
I = (min(detJ, [], 1) ./ max(detJ, [], 1))';
end


function plot_mesh(msh, filename)
figure; clf;
dgmeshplot_curved(msh, 4, 1, 0);
set(findobj(gca, 'Type', 'line'), 'MarkerSize', 3);
set(gcf, 'Position', [114 1 560 560]);
exportgraphics(gcf, filename, 'Resolution', 200);
end


function plot_mesh_hl(msh, red_tris, filename)
figure; clf;
hh = dgmeshplot_curved(msh, 4, 1, 0);
if ~isempty(red_tris)
    set(hh(red_tris), 'FaceColor', [1 0.4 0.4]);
end
set(findobj(gca, 'Type', 'line'), 'MarkerSize', 3);
set(gcf, 'Position', [114 1 560 560]);
exportgraphics(gcf, filename, 'Resolution', 200);
end


function plot_mesh_verts(msh, red_verts, filename)
figure; clf;
dgmeshplot_curved(msh, 4, 1, 0);
set(findobj(gca, 'Type', 'line'), 'MarkerSize', 3);
if ~isempty(red_verts)
    hold on;
    plot(msh.p(1, red_verts), msh.p(2, red_verts), 'o', ...
        'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0], ...
        'MarkerSize', 10);
    hold off;
end
set(gcf, 'Position', [114 1 560 560]);
exportgraphics(gcf, filename, 'Resolution', 200);
end


function m = boundary_linear_mask(msh)
np = size(msh.p, 2);
t = double(msh.t) + 1;
t2t = double(msh.t2t);
m = false(np, 1);
for it = 1:size(t, 2)
    for j = 1:size(t2t, 1)
        if t2t(j, it) < 0
            v = t(:, it); v(j) = [];
            m(v) = true;
        end
    end
end
end


function [ok, pair] = is_3_connected_diagnostic(msh)
% Return [true, []] if graph of msh.t is 3-connected; otherwise return
% [false, [i,j]] a witness pair whose removal disconnects it.
np = size(msh.p, 2);
t  = double(msh.t) + 1;
I = reshape(t([1 2 3 1 2 3], :), [], 1);
J = reshape(t([2 3 1 3 1 2], :), [], 1);
A = sparse(I, J, 1, np, np);
A = (A + A') > 0;
for i = 1:np-1
    for j = i+1:np
        keep = true(np, 1); keep(i) = false; keep(j) = false;
        if ~is_connected_sub(A, keep)
            ok = false; pair = [i, j]; return;
        end
    end
end
ok = true; pair = [];
end


function c = is_connected_sub(A, keep)
verts = find(keep);
if isempty(verts), c = true; return; end
visited = false(size(A, 1), 1);
queue = verts(1);
visited(queue) = true;
while ~isempty(queue)
    v = queue(end); queue(end) = [];
    nbrs = find(A(v, :) & keep');
    for nb = nbrs
        if ~visited(nb)
            visited(nb) = true;
            queue(end+1) = nb; %#ok<AGROW>
        end
    end
end
c = all(visited(verts));
end


function m = boundary_dg_mask(msh)
data = dginit(msh);
[~, edg] = getbndnodes(msh, data);          % linear index into (ns, 1, nt)
[ns, d, nt] = size(msh.p1);
m2 = false(ns, nt);
m2(edg) = true;
m = repmat(reshape(m2, ns, 1, nt), 1, d, 1);
end
