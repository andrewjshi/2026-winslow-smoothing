% stadium_preview.m - Stadium domain: a rectangle of length 2a with a
% half-disk of radius b capping each end. Convex, smoothly curved on the
% half-disk portions, flat on the long sides.
%
% Signed distance function: fd(p) = sqrt(max(|x|-a, 0)^2 + y^2) - b.

rng(3);
porder = 4;
a = 1.0;        % half the straight rectangle length
b = 0.5;        % radius of the capping half-disks (also half-height)

n_bnd = 80;
n_int = 100;

% --- Boundary: sample uniformly by arc length around the stadium ---
perim_straight = 2*a;
perim_semi     = pi*b;
P = 2*perim_straight + 2*perim_semi;
s_vals = P * (0:n_bnd-1)' / n_bnd;
bnd = zeros(n_bnd, 2);
for i = 1:n_bnd
    s = s_vals(i);
    if s < perim_straight
        % bottom edge, from (-a, -b) rightward
        bnd(i, :) = [-a + s, -b];
    elseif s < perim_straight + perim_semi
        % right semicircle, from (a, -b) CCW to (a, b)
        t = (s - perim_straight)/b - pi/2;
        bnd(i, :) = [a + b*cos(t), b*sin(t)];
    elseif s < 2*perim_straight + perim_semi
        % top edge, from (a, b) leftward to (-a, b)
        bnd(i, :) = [a - (s - perim_straight - perim_semi), b];
    else
        % left semicircle, from (-a, b) CCW to (-a, -b)
        t = (s - 2*perim_straight - perim_semi)/b + pi/2;
        bnd(i, :) = [-a + b*cos(t), b*sin(t)];
    end
end

% --- Signed distance function and interior rejection sampling -----------
stadium_fd = @(p) sqrt(max(abs(p(:,1)) - a, 0).^2 + p(:,2).^2) - b;

x_min = -(a + b) - 0.02;  x_max = (a + b) + 0.02;
y_min = -b - 0.02;         y_max = b + 0.02;
int_pts = zeros(0, 2);
while size(int_pts, 1) < n_int
    pt = [x_min + (x_max - x_min) * rand(), ...
          y_min + (y_max - y_min) * rand()];
    if stadium_fd(pt) < -0.02      % safely inside
        int_pts(end+1, :) = pt;    %#ok<AGROW>
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

figdir = 'figures/stadium_preview';
if ~exist(figdir, 'dir'), mkdir(figdir); end

% --- Tutte embedding onto unit disk ---
msh_t = tutte_embedding(msh, 'TargetShape', 'circle');
msh_t = nodealloc(msh_t, porder);

% --- Identify the boundary-sliver clump on the LEFT of the Tutte disk
%     so we can trace it back to its physical origin. ---
tt = double(msh_t.t) + 1;
nt = size(tt, 2);
centroids = zeros(2, nt);
for it = 1:nt
    centroids(:, it) = mean(msh_t.p(:, tt(:, it)), 2);
end
cx = centroids(1, :);  cy = centroids(2, :);
ang = atan2(cy, cx);  rad = sqrt(cx.^2 + cy.^2);
% Left of disk: moderately wide angular band near the boundary ring.
cand = find(abs(ang) >= 3*pi/4 & rad >= 0.85);

% Largest connected component via dual-graph BFS
t2t = double(msh_t.t2t);
is_cand = false(nt, 1);  is_cand(cand) = true;
visited = false(nt, 1);
best_comp = [];
for seed = cand(:)'
    if visited(seed), continue; end
    comp = []; q = seed; visited(seed) = true;
    while ~isempty(q)
        v = q(end); q(end) = [];
        comp(end+1) = v;  %#ok<AGROW>
        for j = 1:size(t2t, 1)
            nb = t2t(j, v) + 1;
            if nb >= 1 && is_cand(nb) && ~visited(nb)
                visited(nb) = true;
                q(end+1) = nb;  %#ok<AGROW>
            end
        end
    end
    if numel(comp) > numel(best_comp), best_comp = comp; end
end
red_tris = best_comp(:);
fprintf('  left clump (Tutte disk): %d triangles\n', numel(red_tris));

plot_hl(msh,   red_tris, fullfile(figdir, 'mesh_clean.png'),        [900 400]);
plot_hl(msh_t, red_tris, fullfile(figdir, 'mesh_tutte_circle.png'), [600 600]);

% --- Shape optimisation (pinned boundary) ---
t0 = tic;
[msh_o, info] = optimize_shape(msh_t);
t_shape_opt = toc(t0);
msh_o = nodealloc(msh_o, porder);
plot_hl(msh_o, red_tris, fullfile(figdir, 'mesh_optimized_circle.png'), [600 600]);
fprintf('  shape-opt: F_init = %.4e, F_final = %.4e, iters = %d, wall = %.2fs\n', ...
    info.F_initial, info.F_final, info.iterations, t_shape_opt);

% --- Winslow hand-off: shape-opt output as reference, stadium boundary
%     as Dirichlet.
curvep1 = msh_o.p1;
is_bnd_dg = boundary_dg_mask_local(msh_o);
curvep1(is_bnd_dg) = msh.p1(is_bnd_dg);
t0 = tic;
msh_w = elliptic_smoothing(msh_o, curvep1, false);
t_winslow_so = toc(t0);
plot_hl(msh_w, red_tris, fullfile(figdir, 'mesh_winslow_circle.png'), [900 400]);
fprintf('  Winslow (after shape-opt): wall = %.2fs\n', t_winslow_so);

% --- Ablation: Winslow directly from raw Tutte (no shape-opt) ---
curvep1_t = msh_t.p1;
is_bnd_dg_t = boundary_dg_mask_local(msh_t);
curvep1_t(is_bnd_dg_t) = msh.p1(is_bnd_dg_t);
try
    t0 = tic;
    msh_wt = elliptic_smoothing(msh_t, curvep1_t, false);
    t_winslow_raw = toc(t0);
    plot_hl(msh_wt, red_tris, fullfile(figdir, 'mesh_winslow_from_tutte.png'), [900 400]);
    fprintf('  Winslow (from raw Tutte): wall = %.2fs\n', t_winslow_raw);
catch err
    fprintf('  Winslow-from-Tutte FAILED: %s\n', err.message);
end

fprintf('stadium preview: nelem = %d, np = %d\n', ...
    size(msh.t, 2), size(msh.p, 2));

% Per-element Liu-Joe eta on the two Winslow outputs.
fprintf('\nPer-element quality on physical meshes:\n');
report_quality('Winslow from shape-opt (panel e)', msh_w);
if exist('msh_wt', 'var')
    report_quality('Winslow from raw Tutte   (panel c)', msh_wt);
end

fprintf('\nTiming summary:\n');
fprintf('  shape-opt           : %.2fs\n', t_shape_opt);
fprintf('  Winslow (shape-opt) : %.2fs\n', t_winslow_so);
if exist('t_winslow_raw', 'var')
    fprintf('  Winslow (raw Tutte) : %.2fs\n', t_winslow_raw);
    fprintf('  shape-opt + Winslow vs Winslow alone: %.2fs vs %.2fs\n', ...
        t_shape_opt + t_winslow_so, t_winslow_raw);
end

fprintf('\nSaved figures to %s/\n', figdir);


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


function report_quality(label, msh)
eta = tri_eta_linear(msh.p, double(msh.t) + 1);
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


function plot_hl(msh, red_tris, filename, figsize)
figure; clf;
hh = dgmeshplot_curved(msh, 4, 1, 0);
set(hh(red_tris), 'FaceColor', [1 0.4 0.4]);
set(findobj(gca, 'Type', 'line'), 'MarkerSize', 1);
axis equal; axis off;
set(gcf, 'Position', [100 100 figsize(1) figsize(2)]);
exportgraphics(gcf, filename, 'Resolution', 200);
close;
end
