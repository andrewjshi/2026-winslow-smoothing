% comb_slide.m - Comb domain through Tutte + sliding shape-opt + beta-
% homotopy Winslow.
%
% Strategy: pin only the four outer corners of the comb's bounding box on
% the unit disk, let every other boundary vertex slide along the disk
% during shape-opt (configuration B in the notes). This lets the optimizer
% redistribute the four tooth-wedges so their tips are not as squeezed as
% the proportional-arc Tutte gives them. Then a proportional-arc pullback
% gives a new physical comb boundary, and beta-homotopy Winslow runs with
% that boundary.

rng(1);
porder = 4;

% --- Comb parameters (must match comb.m) --------------------------------
W       = 0.92;
H_base  = 0.20;
H_top   = 0.70;
margin  = 0.10;
tooth_w = 0.12;
gap_w   = 0.08;
n_teeth = 4;
n_per_edge_density = 30;
n_interior = 400;

% --- Corners (CCW) ------------------------------------------------------
C = [0, 0; W, 0; W, H_base];
x = W - margin;
for k = 1:n_teeth
    C = [C; x, H_base; x, H_top; x - tooth_w, H_top; x - tooth_w, H_base]; %#ok<AGROW>
    x = x - tooth_w - gap_w;
end
C = [C; 0, H_base];

% --- Boundary samples ---------------------------------------------------
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

% --- Interior + 3-connectedness ------------------------------------------
int_pts = zeros(0, 2);
while size(int_pts, 1) < n_interior
    pt = [W * rand(), H_top * rand()];
    if inpolygon(pt(1), pt(2), bnd_pts(:,1), bnd_pts(:,2))
        int_pts = [int_pts; pt]; %#ok<AGROW>
    end
end
[int_pts, anchor_info] = make_3_connected(bnd_pts, int_pts, ...
    'MinSpacing', 0.01, 'NudgeStep', 0.02);
fprintf('make_3_connected: %d iterations, %d anchors added, %d cuts remaining\n', ...
    anchor_info.iters, anchor_info.added, anchor_info.final_cuts);

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

figdir = 'figures/comb_slide';
if ~exist(figdir, 'dir'), mkdir(figdir); end
plot_mesh(msh_clean, fullfile(figdir, 'mesh_clean.png'));
fprintf('clean comb mesh:\n'); print_quality(msh_clean);

% --- Identify outer-rectangle corners as pin vertices on the disk -------
outer_corners = [0, 0; W, 0; W, H_base; 0, H_base];
pin_v = zeros(size(outer_corners, 1), 1);
for k = 1:numel(pin_v)
    [d, idx] = min(vecnorm(msh_clean.p' - outer_corners(k,:), 2, 2));
    if d > 1e-6
        warning('comb_slide: outer corner %d not found within tolerance (d=%.3g)', k, d);
    end
    pin_v(k) = idx;
end
fprintf('pin vertices: %s\n', mat2str(pin_v(:)'));
fprintf('pin positions: \n'); disp(msh_clean.p(:, pin_v)');

% --- Tutte to disk ------------------------------------------------------
msh_t = tutte_embedding(msh_clean, 'TargetShape', 'circle');
msh_t = nodealloc(msh_t, porder);
fprintf('after Tutte:\n'); print_quality(msh_t);
plot_mesh(msh_t, fullfile(figdir, 'mesh_tutte_circle.png'));

% --- Slide-boundary shape-opt -------------------------------------------
[msh_o, info_o] = optimize_shape(msh_t, ...
    'SlideBoundary', 'circle', ...
    'PinBoundaryNodes', pin_v, ...
    'MaxIterations', 500);
msh_o = nodealloc(msh_o, porder);
fprintf('after slide shape-opt:\n'); print_quality(msh_o);
fprintf('  F_init = %.4e, F_final = %.4e, iterations = %d\n', ...
    info_o.F_initial, info_o.F_final, info_o.iterations);
plot_mesh(msh_o, fullfile(figdir, 'mesh_optimized_circle_slide.png'));

% --- Pullback to derive new comb boundary -------------------------------
msh_phys_new = derive_physical_from_slid_circle(msh_o, msh_clean, porder, pin_v(1));
plot_mesh(msh_phys_new, fullfile(figdir, 'mesh_physical_target.png'));

% --- Beta-homotopy Winslow ----------------------------------------------
n_steps = 10;
ws_opts = struct('alpha', 0.5, 'maxiter', 100, 'tol', 1e-5);
try
    [msh_w, hinfo] = winslow_homotopy(msh_o, msh_phys_new, n_steps, ws_opts);
    fprintf('homotopy: %d/%d stages converged. completed=%d.\n', ...
        hinfo.converged_steps, n_steps, hinfo.completed);
    fprintf('after Winslow:\n'); print_quality(msh_w);
    plot_mesh(msh_w, fullfile(figdir, 'mesh_winslow_circle.png'));
catch ME
    fprintf('Winslow FAILED: %s\n', ME.message);
end

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
