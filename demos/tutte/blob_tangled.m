% blob_tangled.m - Tangled "blob" input through the Tutte pipeline,
% tested against both unit-square and unit-circle computational domains.
%
% The physical domain is a smooth non-convex "blob" with radius
%   r(theta) = 1 + 0.3*sin(3*theta) + 0.15*cos(5*theta) + 0.1*sin(7*theta)
% translated / scaled to fit in [-1, 1]^2.

rng(0);
porder = 4;
n_bnd  = 60;    % boundary samples
n_int  = 100;   % interior samples (before rejection)

% --- 1. Build the blob boundary and point set -----------------------------
theta   = 2*pi*(0:n_bnd-1)'/n_bnd;
r       = 1 + 0.30*sin(3*theta) + 0.15*cos(5*theta) + 0.10*sin(7*theta);
bnd_pts = [r.*cos(theta), r.*sin(theta)];

% Centre and scale to fit snugly in [-0.9, 0.9]^2.
bnd_pts = bnd_pts - mean(bnd_pts, 1);
bnd_pts = bnd_pts / max(abs(bnd_pts(:))) * 0.9;

% Interior points: rejection sampling inside the blob polygon.
int_pts = zeros(0, 2);
while size(int_pts, 1) < n_int
    pt = -0.9 + 1.8*rand(1, 2);
    if inpolygon(pt(1), pt(2), bnd_pts(:,1), bnd_pts(:,2))
        int_pts = [int_pts; pt]; %#ok<AGROW>
    end
end

all_pts = [bnd_pts; int_pts];

% Constrained Delaunay respecting the blob boundary.
n_b = size(bnd_pts, 1);
ce  = [(1:n_b)', [(2:n_b)'; 1]];
dt  = delaunayTriangulation(all_pts, ce);
in  = isInterior(dt);
tri = dt.ConnectivityList(in, :);

bndexpr = {'true'};   % one boundary label for everything
msh_clean = ml2msh(all_pts, tri, bndexpr);
msh_clean = mshcurved(msh_clean, []);
msh_clean = nodealloc(msh_clean, porder);

figdir = 'figures/blob_tangled';
if ~exist(figdir, 'dir'), mkdir(figdir); end

fprintf('clean blob mesh:\n'); print_quality(msh_clean);
[ok, pair] = is_3_connected_diagnostic(msh_clean);
if ok
    fprintf('graph is 3-connected.\n');
else
    fprintf('graph NOT 3-connected; remove %d, %d disconnects at (%.3f, %.3f), (%.3f, %.3f)\n', ...
            pair(1), pair(2), msh_clean.p(1, pair(1)), msh_clean.p(2, pair(1)), ...
            msh_clean.p(1, pair(2)), msh_clean.p(2, pair(2)));
end
plot_mesh(msh_clean, fullfile(figdir, 'mesh_clean.png'));

% --- 2. Tangle violently ----------------------------------------------------
msh_tangled = msh_clean;

is_bnd_lin = boundary_linear_mask(msh_tangled);
amp = 0.4;
d_lin = amp * (rand(size(msh_tangled.p)) - 0.5) * 2;
d_lin(:, is_bnd_lin) = 0;
msh_tangled.p = msh_tangled.p + d_lin;
msh_tangled = nodealloc(msh_tangled, porder);

is_bnd_dg = boundary_dg_mask(msh_tangled);
d_ho = 0.12 * (rand(size(msh_tangled.p1)) - 0.5) * 2;
d_ho(is_bnd_dg) = 0;
msh_tangled.p1 = msh_tangled.p1 + d_ho;

fprintf('\ntangled input:\n'); print_quality(msh_tangled);
plot_mesh(msh_tangled, fullfile(figdir, 'mesh_tangled.png'));

% --- 3. Pipeline with two targets -------------------------------------------
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

    % Winslow: optimised Tutte as reference, tangled input's blob
    % boundary as Dirichlet data.
    curvep1 = msh_o.p1;
    is_bnd_dg_ref = boundary_dg_mask(msh_o);
    curvep1(is_bnd_dg_ref) = msh_tangled.p1(is_bnd_dg_ref);
    msh_w = elliptic_smoothing(msh_o, curvep1, false);

    fprintf('after Winslow:\n'); print_quality(msh_w);
    plot_mesh(msh_w, fullfile(figdir, sprintf('mesh_winslow_%s.png', ts)));
end

fprintf('\nSaved figures to %s/\n', figdir);


% ===========================================================================
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
set(gcf, 'Position', [114 1 560 560]);
exportgraphics(gcf, filename, 'Resolution', 200);
end


function [ok, pair] = is_3_connected_diagnostic(msh)
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
