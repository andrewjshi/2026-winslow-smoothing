% combdome.m -- Comb domain through Tutte to the "dome" reference shape:
% rectangle [0,W] x [0,H_base] capped by a half-ellipse with semi-axes
% (W/2, b), where b is chosen so the cap arc length equals the comb's
% teeth-region arc length. The comb teeth get ~50 degrees of arc each
% on the cap instead of being squeezed into a small disk arc.
%
% Stage 1: just Tutte. No shape opt, no Winslow yet.

rng(1);
porder = 4;

% --- Comb parameters (must match comb.m) --------------------------------
W       = 0.92;
H_base  = 0.20;
H_top   = 0.38;
margin  = 0.10;
tooth_w = 0.18;
gap_w   = 0.10;
n_teeth = 3;
n_per_edge_density = 25;
n_interior = 350;

% --- Comb corners (CCW) -------------------------------------------------
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

% --- Boundary samples with per-point bookkeeping ------------------------
bnd_pts = [];
edge_id = [];   % which corner-edge each boundary sample lives on
edge_t  = [];   % parameter t in [0,1) along that edge
for i = 1:nc
    L = edge_len(i);
    npts = max(1, round(L * n_per_edge_density));
    for k = 0:npts-1
        t = k / npts;
        bnd_pts = [bnd_pts; C(i, :) + t * edge_vec(i, :)]; %#ok<AGROW>
        edge_id = [edge_id; i]; %#ok<AGROW>
        edge_t  = [edge_t;  t]; %#ok<AGROW>
    end
end
n_b = size(bnd_pts, 1);

% --- Interior samples + 3-connectedness ---------------------------------
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
ce  = [(1:n_b)', [(2:n_b)'; 1]];
dt  = delaunayTriangulation(all_pts, ce);
in  = isInterior(dt);
tri = dt.ConnectivityList(in, :);

bndexpr = {'true'};
msh_clean = ml2msh(all_pts, tri, bndexpr);
msh_clean = mshcurved(msh_clean, []);
msh_clean = nodealloc(msh_clean, porder);

figdir = 'figures/combdome';
if ~exist(figdir, 'dir'), mkdir(figdir); end
plot_mesh(msh_clean, fullfile(figdir, 'mesh_clean.png'));
fprintf('clean comb mesh:\n'); print_quality(msh_clean);

% --- Dome geometry ------------------------------------------------------
L_teeth = sum(edge_len(3:nc-1));
a_dome  = W/2;
b_dome  = solve_b_for_arc(L_teeth, a_dome);
fprintf('dome: a=%.4f, b=%.4f, half-perim=%.4f (target %.4f)\n', ...
    a_dome, b_dome, ellipse_half_perim(a_dome, b_dome), L_teeth);

% Precompute (theta, s(theta)) on the cap for arc-length-preserving mapping.
N_theta    = 4000;
theta_grid = linspace(0, pi, N_theta+1)';
ds_dtheta  = sqrt(a_dome^2 * sin(theta_grid).^2 + b_dome^2 * cos(theta_grid).^2);
s_grid     = [0; cumsum((ds_dtheta(1:end-1) + ds_dtheta(2:end))/2 * (pi/N_theta))];
s_grid     = s_grid * (L_teeth / s_grid(end));   % renormalize to exact L_teeth

% --- Map every boundary mesh vertex onto the dome -----------------------
target_bnd = zeros(n_b, 2);
for k = 1:n_b
    i = edge_id(k);
    if i == 1 || i == 2 || i == nc
        % Edges 1 (bottom), 2 (right), nc (left) are shared with the dome
        % rectangle exactly, so the boundary points stay put.
        target_bnd(k, :) = bnd_pts(k, :);
    else
        % Teeth-region edge: place on the cap by arc length from corner 3.
        s = sum(edge_len(3:i-1)) + edge_t(k) * edge_len(i);
        theta = interp1(s_grid, theta_grid, s, 'linear');
        target_bnd(k, 1) = W/2 + a_dome * cos(theta);
        target_bnd(k, 2) = H_base + b_dome * sin(theta);
    end
end

% --- Tutte solve (preserve the dome boundary we just placed) ------------
msh_pre = msh_clean;
msh_pre.p(1, 1:n_b) = target_bnd(:, 1)';
msh_pre.p(2, 1:n_b) = target_bnd(:, 2)';

msh_t = tutte_embedding(msh_pre, 'PreserveBoundary', true);
msh_t = nodealloc(msh_t, porder);
fprintf('after Tutte (dome reference):\n'); print_quality(msh_t);
plot_mesh(msh_t, fullfile(figdir, 'mesh_tutte_dome.png'));

% --- Shape-opt with pinned dome boundary --------------------------------
[msh_o, info_o] = optimize_shape(msh_t);
msh_o = nodealloc(msh_o, porder);
fprintf('after shape-opt:\n'); print_quality(msh_o);
fprintf('  F_init = %.4e, F_final = %.4e, iterations = %d\n', ...
    info_o.F_initial, info_o.F_final, info_o.iterations);
plot_mesh(msh_o, fullfile(figdir, 'mesh_optimized_dome.png'));

% --- Winslow (no homotopy): dome -> comb in one shot --------------------
% The dome reference is close enough to the comb that the elliptic solver
% converges directly. Beta-homotopy interpolates DG boundary nodes
% linearly between dome and comb, which cuts across re-entrant corners
% at intermediate beta and locks in inverted curved elements there.
ws_opts = struct('alpha', 0.5, 'maxiter', 100, 'tol', 1e-5);

fprintf('\n=== Winslow (no homotopy): raw Tutte -> comb ===\n');
try
    [msh_w_t, ~] = winslow_homotopy(msh_t, msh_clean, 1, ws_opts);
    fprintf('after Winslow (from raw Tutte):\n'); print_quality(msh_w_t);
    plot_mesh(msh_w_t, fullfile(figdir, 'mesh_winslow_from_tutte.png'));
catch ME
    fprintf('Winslow from raw Tutte FAILED: %s\n', ME.message);
end

fprintf('\n=== Winslow (no homotopy): optimized Tutte -> comb ===\n');
try
    [msh_w_o, ~] = winslow_homotopy(msh_o, msh_clean, 1, ws_opts);
    fprintf('after Winslow (from optimized Tutte):\n'); print_quality(msh_w_o);
    plot_mesh(msh_w_o, fullfile(figdir, 'mesh_winslow_from_optimized.png'));
catch ME
    fprintf('Winslow from optimized Tutte FAILED: %s\n', ME.message);
end

fprintf('\nSaved figures to %s/\n', figdir);


% =========================================================================
function p = ellipse_half_perim(a, b)
p = pi * (3*(a+b) - sqrt((3*a+b)*(a+3*b))) / 2;
end

function b = solve_b_for_arc(L, a)
lo = 1e-6;
hi = max(1, 4*L/pi);
while ellipse_half_perim(a, hi) < L
    hi = 2*hi;
end
for it = 1:80
    mid = (lo + hi)/2;
    if ellipse_half_perim(a, mid) < L
        lo = mid;
    else
        hi = mid;
    end
end
b = (lo + hi)/2;
end

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
