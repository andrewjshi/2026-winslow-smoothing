% starfish_preview.m - Stylized five-armed starfish boundary + constrained
% Delaunay mesh. Preview only; pipeline hand-off is added once the shape
% is approved.
%
% Parametric form:
%   r(theta) = 1 + 0.40 * cos(5*theta) + 0.03 * sin(2*theta)
%   x = r * cos(theta), y = r * sin(theta)
% The principal cos(5*theta) term gives 5 arms; the small sin(2*theta)
% breaks exact 5-fold symmetry so the figure looks organic rather than
% CAD-ish.

rng(7);
porder = 4;
n_bnd  = 60;
n_int  = 80;

t = 2*pi * (0:n_bnd-1)' / n_bnd;
r = 1 + 0.40 * cos(5*t) + 0.03 * sin(2*t);
xb = r .* cos(t);
yb = r .* sin(t);
bnd = [xb, yb];

% Signed distance function for the starfish boundary (radial
% approximation: accurate near the boundary, sufficient for the DG-node
% projection that mshcurved/nodealloc perform).
starfish_fd = @(p) sqrt(sum(p.^2, 2)) - ( ...
    1 + 0.40 * cos(5 * atan2(p(:,2), p(:,1))) ...
      + 0.03 * sin(2 * atan2(p(:,2), p(:,1))));

% Interior rejection sampling
x_min = min(xb) - 0.05; x_max = max(xb) + 0.05;
y_min = min(yb) - 0.05; y_max = max(yb) + 0.05;
int_pts = zeros(0, 2);
while size(int_pts, 1) < n_int
    pt = [x_min + (x_max - x_min) * rand(), ...
          y_min + (y_max - y_min) * rand()];
    if inpolygon(pt(1), pt(2), xb, yb)
        int_pts(end+1, :) = pt;  %#ok<AGROW>
    end
end

all_pts = [bnd; int_pts];
n_b = size(bnd, 1);
ce = [(1:n_b)', [(2:n_b)'; 1]];
dt = delaunayTriangulation(all_pts, ce);
in = isInterior(dt);
tri = dt.ConnectivityList(in, :);

bndexpr = {'true'};
msh = ml2msh(all_pts, tri, bndexpr, starfish_fd, {});
msh = mshcurved(msh, 1);
msh = nodealloc(msh, porder);

figdir = 'figures/starfish_preview';
if ~exist(figdir, 'dir'), mkdir(figdir); end

% --- Tutte embedding on unit circle target ---
msh_t = tutte_embedding(msh, 'TargetShape', 'circle');
msh_t = nodealloc(msh_t, porder);

% --- Shape optimisation (pinned boundary) ---
[msh_o, info] = optimize_shape(msh_t);
msh_o = nodealloc(msh_o, porder);

% --- Select the bottom-left arm of the PHYSICAL starfish. This arm
%     has the worst-shaped triangles (long skinny slivers along the
%     boundary), making it the natural "hard case" to track. Since the
%     topology is unchanged through Tutte + shape-opt + Winslow, the
%     same triangle indices can be highlighted in every stage figure.
t = double(msh.t) + 1;
nt = size(t, 2);
centroids = zeros(2, nt);
for it = 1:nt
    centroids(:, it) = mean(msh.p(:, t(:, it)), 2);
end
cx = centroids(1, :);  cy = centroids(2, :);
ang = atan2(cy, cx);   rad = sqrt(cx.^2 + cy.^2);
% Bottom-left arm: centroid angle roughly in [-pi, -3*pi/5] (= -180 to
% -108 deg), past the body-centre radius (lowered to 0.6 to include the
% base triangles that otherwise leave a concave notch in the highlight).
cand = find(ang >= -pi & ang <= -3*pi/5 & rad >= 0.6);
% Filter to the largest connected component via dual-graph BFS.
t2t = double(msh_t.t2t);
is_cand = false(nt, 1);  is_cand(cand) = true;
visited = false(nt, 1);
best_comp = [];
for seed = cand(:)'
    if visited(seed), continue; end
    % BFS over the subgraph restricted to cand triangles
    comp = []; queue = seed; visited(seed) = true;
    while ~isempty(queue)
        v = queue(end); queue(end) = [];
        comp(end+1) = v;  %#ok<AGROW>
        for j = 1:size(t2t, 1)
            nb = t2t(j, v) + 1;
            if nb >= 1 && is_cand(nb) && ~visited(nb)
                visited(nb) = true;
                queue(end+1) = nb;  %#ok<AGROW>
            end
        end
    end
    if numel(comp) > numel(best_comp), best_comp = comp; end
end
red_tris = best_comp(:);
fprintf('  red arm (physical bottom-left): %d triangles\n', numel(red_tris));

% (Single-arm plots removed; the two-arm versions below are the only
% figures produced for this workflow.)

% --- Winslow hand-off: shape-opt output is reference, physical starfish
%     boundary is Dirichlet data.
curvep1 = msh_o.p1;
is_bnd_dg = boundary_dg_mask_local(msh_o);
curvep1(is_bnd_dg) = msh.p1(is_bnd_dg);
msh_w = elliptic_smoothing(msh_o, curvep1, false);


% --- Ablation: Winslow from raw Tutte output (no shape-opt). Is the
%     shape-opt intermediate actually contributing, or would Winslow
%     alone be enough?
curvep1_t = msh_t.p1;
is_bnd_dg_t = boundary_dg_mask_local(msh_t);
curvep1_t(is_bnd_dg_t) = msh.p1(is_bnd_dg_t);
try
    msh_wt = elliptic_smoothing(msh_t, curvep1_t, false);
    fprintf('  Winslow-from-Tutte: converged.\n');
catch err
    fprintf('  Winslow-from-Tutte FAILED: %s\n', err.message);
end

% Plot each pipeline stage with the red (bottom-left) arm highlighted.
plot_hl(msh,   red_tris, fullfile(figdir, 'mesh_clean.png'));
plot_hl(msh_t, red_tris, fullfile(figdir, 'mesh_tutte_circle.png'));
plot_hl(msh_o, red_tris, fullfile(figdir, 'mesh_optimized_circle.png'));
plot_hl(msh_w, red_tris, fullfile(figdir, 'mesh_winslow_circle.png'));
if exist('msh_wt', 'var')
    plot_hl(msh_wt, red_tris, fullfile(figdir, 'mesh_winslow_from_tutte.png'));
end

fprintf('starfish preview:\n');
fprintf('  clean: nelem = %d, np = %d\n', size(msh.t, 2), size(msh.p, 2));
fprintf('  shape-opt: F_init = %.4e, F_final = %.4e, iters = %d\n', ...
    info.F_initial, info.F_final, info.iterations);

% Compare per-element mean eta on the physical meshes we produced.
fprintf('\nPer-element Liu-Joe eta on physical meshes:\n');
report_quality('clean input             ', msh);
if exist('msh_w', 'var')
    report_quality('pipeline (with shape-opt)', msh_w);
end
if exist('msh_wt', 'var')
    report_quality('pipeline (raw Tutte)    ', msh_wt);
end

fprintf('\nSaved figures to %s/\n', figdir);


% =========================================================================
function eta = tri_eta_linear(p, t)
% Liu-Joe eta on each linear triangle (msh.p, msh.t one-based).
nt = size(t, 2);
eta = zeros(nt, 1);
A1 = 4/3; A2 = -2/3; A3 = 4/3;
det_Minv = 2/sqrt(3);
for it = 1:nt
    p1 = p(:, t(1, it));  p2 = p(:, t(2, it));  p3 = p(:, t(3, it));
    a = p2(1) - p1(1);  b = p3(1) - p1(1);
    c = p2(2) - p1(2);  d = p3(2) - p1(2);
    detE = a * d - b * c;
    if detE <= 0
        eta(it) = Inf;
        continue;
    end
    Q = (a^2 + c^2) * A1 + 2 * (a*b + c*d) * A2 + (b^2 + d^2) * A3;
    D = 2 * det_Minv * detE;
    eta(it) = Q / D;
end
end


function report_quality(label, msh)
t = double(msh.t) + 1;
eta = tri_eta_linear(msh.p, t);
fprintf('  %s: mean(eta^2)=%.3f, max(eta)=%.3f, sum(eta^2)=%.1f\n', ...
    label, mean(eta.^2), max(eta), sum(eta.^2));
end


function m = boundary_dg_mask_local(msh)
data = dginit(msh);
[~, edg] = getbndnodes(msh, data);
[ns, d, nt] = size(msh.p1);
m2 = false(ns, nt);
m2(edg) = true;
m = repmat(reshape(m2, ns, 1, nt), 1, d, 1);
end


function plot_hl(msh, red_tris, filename)
figure; clf;
hh = dgmeshplot_curved(msh, 4, 1, 0);
set(hh(red_tris), 'FaceColor', [1 0.4 0.4]);
set(findobj(gca, 'Type', 'line'), 'MarkerSize', 1);
axis equal; axis off;
set(gcf, 'Position', [100 100 600 600]);
exportgraphics(gcf, filename, 'Resolution', 200);
close;
end


function plot_hl2(msh, red_tris, green_tris, filename)
figure; clf;
hh = dgmeshplot_curved(msh, 4, 1, 0);
set(hh(red_tris),   'FaceColor', [1 0.4 0.4]);
set(hh(green_tris), 'FaceColor', [0.4 0.8 0.4]);
set(findobj(gca, 'Type', 'line'), 'MarkerSize', 1);
axis equal; axis off;
set(gcf, 'Position', [100 100 600 600]);
exportgraphics(gcf, filename, 'Resolution', 200);
close;
end
