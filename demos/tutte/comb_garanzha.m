% comb_garanzha.m - Comb pipeline using the Garanzha 2021 foldover-free
% penalty method (Escobar-regularized barrier) in place of Picard Winslow.
%
% Setup is identical to comb.m through the Tutte and shape-opt stages.
% The Winslow stage is replaced by winslow_garanzha which optimizes
% interior linear-vertex positions to drive every triangle's det J to a
% positive value, while keeping the boundary fixed at the comb's physical
% positions.

rng(1);
porder = 4;

% Comb parameters (matches the current comb.m)
W       = 0.92;
H_base  = 0.20;
H_top   = 0.38;
margin  = 0.10;
tooth_w = 0.18;
gap_w   = 0.10;
n_teeth = 3;
n_per_edge_density = 25;
n_interior = 350;

% --- Build comb (same steps as comb.m) ----------------------------------
C = [0, 0; W, 0; W, H_base];
x = W - margin;
for k = 1:n_teeth
    C = [C; x, H_base; x, H_top; x - tooth_w, H_top; x - tooth_w, H_base]; %#ok<AGROW>
    x = x - tooth_w - gap_w;
end
C = [C; 0, H_base];

nc = size(C, 1);
edge_vec = C([2:nc 1], :) - C;
edge_len = vecnorm(edge_vec, 2, 2);
bnd_pts = [];
for i = 1:nc
    L = edge_len(i);
    npts = max(1, round(L * n_per_edge_density));
    for k = 0:npts-1
        t = k / npts;
        bnd_pts = [bnd_pts; C(i, :) + t * edge_vec(i, :)]; %#ok<AGROW>
    end
end

int_pts = zeros(0, 2);
while size(int_pts, 1) < n_interior
    pt = [W * rand(), H_top * rand()];
    if inpolygon(pt(1), pt(2), bnd_pts(:,1), bnd_pts(:,2))
        int_pts = [int_pts; pt]; %#ok<AGROW>
    end
end
[int_pts, anchor_info] = make_3_connected(bnd_pts, int_pts, ...
    'MinSpacing', 0.01, 'NudgeStep', 0.02);
fprintf('make_3_connected: %d iterations, %d anchors added\n', ...
    anchor_info.iters, anchor_info.added);

all_pts = [bnd_pts; int_pts];
n_b = size(bnd_pts, 1);
ce  = [(1:n_b)', [(2:n_b)'; 1]];
dt  = delaunayTriangulation(all_pts, ce);
in  = isInterior(dt);
tri = dt.ConnectivityList(in, :);

bndexpr = {'true'};
msh_clean = ml2msh(all_pts, tri, bndexpr);
msh_clean = mshcurved(msh_clean, []);
msh_clean = nodealloc(msh_clean, porder);

figdir = 'figures/comb_garanzha';
if ~exist(figdir, 'dir'), mkdir(figdir); end
plot_mesh(msh_clean, fullfile(figdir, 'mesh_clean.png'));
fprintf('clean comb mesh:\n'); print_quality(msh_clean);

% --- Tutte to disk ------------------------------------------------------
msh_t = tutte_embedding(msh_clean, 'TargetShape', 'circle');
msh_t = nodealloc(msh_t, porder);
fprintf('after Tutte:\n'); print_quality(msh_t);
plot_mesh(msh_t, fullfile(figdir, 'mesh_tutte_circle.png'));

% --- Shape-opt on disk --------------------------------------------------
[msh_o, info] = optimize_shape(msh_t);
msh_o = nodealloc(msh_o, porder);
fprintf('after shape-opt:\n'); print_quality(msh_o);
fprintf('  F_init = %.4e, F_final = %.4e, iters = %d\n', ...
    info.F_initial, info.F_final, info.iterations);
plot_mesh(msh_o, fullfile(figdir, 'mesh_optimized_circle.png'));

% --- Garanzha foldover-free smoothing in place of Winslow ---------------
% Reference: shape-opt'd disk mesh's linear-vertex positions.
% Target boundary: linear-vertex positions on the physical comb.
is_bnd = boundary_mask_local(msh_o);
bnd_v  = find(is_bnd);
bnd_target = msh_clean.p(:, bnd_v);

[msh_w_lin, ginfo] = winslow_garanzha(msh_o, bnd_target, ...
    'Lambda',   1.0, ...
    'MaxOuter', 30, ...
    'MaxInner', 200, ...
    'EpsInit',  1.0, ...
    'Display',  'iter');

fprintf('garanzha: converged=%d.\n', ginfo.converged);
fprintf('  outer-iter min(detJ) trace: ');
fprintf('%+.2e ', ginfo.min_detJ_per_outer(~isnan(ginfo.min_detJ_per_outer)));
fprintf('\n  outer-iter energy trace:    ');
fprintf('%.3e ', ginfo.energy_per_outer(~isnan(ginfo.energy_per_outer)));
fprintf('\n');

% Re-allocate high-order DG nodes from updated linear vertices on a
% straight-sided boundary. We start from msh_clean (which has the comb
% topology + comb fd cleared to straight), copy in Garanzha's linear
% vertex positions (interior + boundary unchanged), then nodealloc.
msh_w = msh_clean;
msh_w.p = msh_w_lin.p;
msh_w = nodealloc(msh_w, porder);
fprintf('after garanzha (linear-only):\n'); print_quality(msh_w);
plot_mesh(msh_w, fullfile(figdir, 'mesh_garanzha_circle.png'));

fprintf('\nSaved figures to %s/\n', figdir);


% =========================================================================
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
axis equal; axis off;
set(gcf, 'Position', [114 1 800 500]);
exportgraphics(gcf, filename, 'Resolution', 200);
close;
end


function is_bnd = boundary_mask_local(msh)
np = size(msh.p, 2);
t = double(msh.t) + 1;
t2t = double(msh.t2t);
is_bnd = false(np, 1);
for it = 1:size(t, 2)
    for j = 1:size(t2t, 1)
        if t2t(j, it) < 0
            v = t(:, it); v(j) = [];
            is_bnd(v) = true;
        end
    end
end
end
