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

fprintf('clean L-shape mesh:\n'); print_quality(msh_clean);
plot_mesh(msh_clean, fullfile(figdir, 'mesh_clean.png'));

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
% Step 3: Tutte + shape-opt pipeline, two target shapes.
% ---------------------------------------------------------------------
for target_shape = {'square', 'circle'}
    ts = target_shape{1};
    fprintf('\n=== target: %s ===\n', ts);

    msh_t = tutte_embedding(msh_tangled, 'TargetShape', ts);
    msh_t = nodealloc(msh_t, porder);

    fprintf('after Tutte:\n'); print_quality(msh_t);
    plot_mesh(msh_t, fullfile(figdir, sprintf('mesh_tutte_%s.png', ts)));

    [msh_o, info] = optimize_shape(msh_t);
    msh_o = nodealloc(msh_o, porder);

    fprintf('after shape-opt:\n'); print_quality(msh_o);
    fprintf('  F_init = %.4e, F_final = %.4e, iterations = %d\n', ...
            info.F_initial, info.F_final, info.iterations);
    plot_mesh(msh_o, fullfile(figdir, sprintf('mesh_optimized_%s.png', ts)));

    % Winslow: optimised Tutte as reference, tangled input's L-shape
    % boundary as Dirichlet data.
    curvep1 = msh_o.p1;
    is_bnd_dg_ref = boundary_dg_mask(msh_o);
    curvep1(is_bnd_dg_ref) = msh_tangled.p1(is_bnd_dg_ref);
    msh_w = elliptic_smoothing(msh_o, curvep1, false);

    fprintf('after Winslow:\n'); print_quality(msh_w);
    plot_mesh(msh_w, fullfile(figdir, sprintf('mesh_winslow_%s.png', ts)));
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
