% unstructured_tutte.m - Tutte + shape-optimisation pipeline applied to
% a mildly-tangled unstructured Delaunay triangulation of the unit
% square. Random interior points + uniformly-spaced boundary points
% give varying valences, so Tutte must solve a non-trivial Laplacian
% and the shape optimiser has real work to do.

rng(42);              % reproducibility
porder       = 4;
n_per_side   = 10;    % boundary-point spacing
n_interior   = 60;    % random interior points

% --- 1. Generate points ---------------------------------------------------
% Boundary: uniform spacing around the unit-square perimeter.
seg = (0 : n_per_side - 1)' / n_per_side;        % open intervals
boundary_pts = [[seg, zeros(n_per_side, 1)];                % bottom
                [ones(n_per_side, 1), seg];                  % right
                [1 - seg, ones(n_per_side, 1)];              % top
                [zeros(n_per_side, 1), 1 - seg]];            % left

% Interior: random in (0.05, 0.95)^2 so nothing lands near the boundary.
interior_pts = 0.05 + 0.9 * rand(n_interior, 2);

all_pts = [boundary_pts; interior_pts];

% --- 2. Delaunay triangulate ---------------------------------------------
tri = delaunay(all_pts(:, 1), all_pts(:, 2));

% --- 3. Wrap as a dgmatlab msh struct ------------------------------------
bndexpr = {
    'all(p(:,2) < 1e-6)',         % bottom (y = 0)
    'all(p(:,1) > 1 - 1e-6)',     % right  (x = 1)
    'all(p(:,2) > 1 - 1e-6)',     % top    (y = 1)
    'all(p(:,1) < 1e-6)'          % left   (x = 0)
};
msh_input = ml2msh(all_pts, tri, bndexpr);
msh_input = mshcurved(msh_input, []);
msh_input = nodealloc(msh_input, porder);

% --- Stage snapshots -----------------------------------------------------
figdir = 'figures/unstructured_tutte';
if ~exist(figdir, 'dir'), mkdir(figdir); end

fprintf('clean Delaunay:\n'); print_quality(msh_input);
plot_mesh(msh_input, fullfile(figdir, 'mesh_clean.png'));

% --- 4. Mild tangle of interior nodes ------------------------------------
msh_input_tangled = msh_input;
is_bnd_lin = boundary_linear_mask(msh_input_tangled);
amp = 0.1;
d_lin = amp * (rand(size(msh_input_tangled.p)) - 0.5) * 2;
d_lin(:, is_bnd_lin) = 0;
msh_input_tangled.p = msh_input_tangled.p + d_lin;
msh_input_tangled = nodealloc(msh_input_tangled, porder);

is_bnd_dg = boundary_dg_mask(msh_input_tangled);
d_ho = 0.03 * (rand(size(msh_input_tangled.p1)) - 0.5) * 2;
d_ho(is_bnd_dg) = 0;
msh_input_tangled.p1 = msh_input_tangled.p1 + d_ho;

fprintf('\nmildly tangled input:\n'); print_quality(msh_input_tangled);
plot_mesh(msh_input_tangled, fullfile(figdir, 'mesh_tangled.png'));

% --- 5. Tutte embedding ---------------------------------------------------
% Preserve input boundary positions (they already sit on the unit square,
% and we want the graph neighbours of each boundary vertex to connect to
% the same physical location they did in the input).
msh_tutte = tutte_embedding(msh_input_tangled, 'PreserveBoundary', true);
msh_tutte = nodealloc(msh_tutte, porder);

fprintf('\nafter Tutte:\n'); print_quality(msh_tutte);
plot_mesh(msh_tutte, fullfile(figdir, 'mesh_tutte.png'));

% --- 6. Shape optimisation -----------------------------------------------
[msh_opt, info] = optimize_shape(msh_tutte, 'Display', 'iter');
msh_opt = nodealloc(msh_opt, porder);

fprintf('\nafter shape optimisation:\n'); print_quality(msh_opt);
fprintf('  F_init = %.6e, F_final = %.6e, iterations = %d\n', ...
        info.F_initial, info.F_final, info.iterations);
plot_mesh(msh_opt, fullfile(figdir, 'mesh_optimized.png'));

% --- 7. Winslow smoothing: optimised Tutte as reference,
%        tangled input as starting config (its boundary is Dirichlet) ---
msh_winslow = elliptic_smoothing(msh_opt, msh_input_tangled.p1, false);

fprintf('\nafter Winslow:\n'); print_quality(msh_winslow);
plot_mesh(msh_winslow, fullfile(figdir, 'mesh_winslow.png'));

fprintf('\nSaved figures to %s/\n', figdir);


function print_quality(msh)
[eta, I] = quality(msh);
fprintf('  nelem = %d\n', numel(eta));
fprintf('  eta   : min = %.4g, mean = %.4g, max = %.4g\n', ...
        min(eta), mean(eta), max(eta));
fprintf('  I     : min = %.4g (I<=0 => inverted)\n', min(I));
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
set(gcf, 'Position', [114 1 560 420]);
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


function m = boundary_dg_mask(msh)
data = dginit(msh);
[~, edg] = getbndnodes(msh, data);
[ns, d, nt] = size(msh.p1);
m2 = false(ns, nt);
m2(edg) = true;
m = repmat(reshape(m2, ns, 1, nt), 1, d, 1);
end
