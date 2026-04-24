% annulus_tutte.m - Cut-to-disk + Tutte pipeline applied to an annular
% domain (unit disk with an off-centre hole). Demonstrates that a
% multiply connected input can be slit into a topological disk via
% cut_to_disk.m, after which Tutte's algorithm handles it unchanged.

rng(3);
porder = 4;

% --- 1. Build the annular physical domain -------------------------------
% Outer: circle radius R_out centred at origin. Inner: circle radius
% R_in centred at (0.2, 0.0).
R_out = 1.0;  n_out = 48;
R_in  = 0.30; n_in  = 18;
c_in  = [0.2; 0.0];

theta_o = 2*pi*(0:n_out-1)'/n_out;
outer = [R_out*cos(theta_o), R_out*sin(theta_o)];

theta_i = 2*pi*(0:n_in-1)'/n_in;
% The inner boundary must appear CLOCKWISE for Delaunay to treat the
% inside as a hole.  We reverse the theta sweep.
inner = [c_in(1) + R_in*cos(-theta_i), c_in(2) + R_in*sin(-theta_i)];

% Constraint edges: outer loop closes on itself, inner loop closes on
% itself. Inner indices must come after outer in the point list.
ce_outer = [(1:n_out)', [(2:n_out)'; 1]];
ce_inner = n_out + [(1:n_in)', [(2:n_in)'; 1]];

% Interior points: rejection-sample inside the annulus.
n_int = 120;
int_pts = zeros(0, 2);
while size(int_pts, 1) < n_int
    pt = -R_out + 2*R_out*rand(1, 2);
    if norm(pt) < R_out - 0.05 && norm(pt - c_in') > R_in + 0.05
        int_pts(end+1, :) = pt; %#ok<AGROW>
    end
end

all_pts = [outer; inner; int_pts];
ce = [ce_outer; ce_inner];

dt = delaunayTriangulation(all_pts, ce);
in = isInterior(dt);
tri = dt.ConnectivityList(in, :);

bndexpr = {'true'};
msh_annulus = ml2msh(all_pts, tri, bndexpr);
msh_annulus = mshcurved(msh_annulus, []);

figdir = 'figures/annulus_tutte';
if ~exist(figdir, 'dir'), mkdir(figdir); end

fprintf('annular mesh:\n');
fprintf('  nelem = %d, np = %d\n', size(msh_annulus.t, 2), size(msh_annulus.p, 2));
plot_linear(msh_annulus, fullfile(figdir, 'mesh_annulus.png'));

% --- 2. Cut-to-disk ----------------------------------------------------
[msh_disk, info] = cut_to_disk(msh_annulus);
fprintf('\nafter cut_to_disk:\n');
fprintf('  n_holes detected = %d\n', info.n_holes);
for k = 1:numel(info.cuts)
    fprintf('  cut %d: %d vertices, from v=%d to v=%d\n', ...
            k, numel(info.cuts{k}), info.cuts{k}(1), info.cuts{k}(end));
end
fprintf('  nelem = %d, np = %d (vs %d originally)\n', ...
        size(msh_disk.t, 2), size(msh_disk.p, 2), size(msh_annulus.p, 2));

% Verify: one boundary loop.
boundary_faces = sum(sum(double(msh_disk.t2t) < 0));
fprintf('  boundary faces = %d\n', boundary_faces);

% Explicit loop count
loops = debug_extract_loops(msh_disk);
fprintf('  # boundary loops in msh_disk = %d\n', numel(loops));
for k = 1:numel(loops)
    fprintf('    loop %d: %d vertices\n', k, numel(loops{k}));
end

plot_linear_with_cut(msh_disk, info.cuts{1}, ...
                     fullfile(figdir, 'mesh_slit.png'));

% --- 3. Tutte + shape-opt on the slit mesh -----------------------------
msh_disk_p4 = nodealloc(msh_disk, porder);
plot_linear(msh_disk_p4, fullfile(figdir, 'mesh_slit_p4.png'));

msh_tutte = tutte_embedding(msh_disk_p4, 'TargetShape', 'circle');
msh_tutte = nodealloc(msh_tutte, porder);
fprintf('\nafter Tutte (circle target):\n');
fprintf('  nelem = %d\n', size(msh_tutte.t, 2));
plot_linear(msh_tutte, fullfile(figdir, 'mesh_tutte_circle.png'));

[msh_opt, optinfo] = optimize_shape(msh_tutte);
msh_opt = nodealloc(msh_opt, porder);
fprintf('\nafter shape-opt:\n  iterations = %d\n', optinfo.iterations);
plot_linear(msh_opt, fullfile(figdir, 'mesh_optimized_circle.png'));

% --- 4. Winslow: map the optimized disk reference back to the annulus -----
% Reference mesh: msh_opt (topology = slit disk, positions = unit disk).
% Dirichlet data: boundary DG nodes should sit at the original annular
% positions, which are exactly what msh_disk_p4.p1 holds (cut_to_disk
% preserved physical positions; only the unit-disk mapping from Tutte
% moved them).
curvep1 = msh_opt.p1;
is_bnd_dg = boundary_dg_mask(msh_opt);
curvep1(is_bnd_dg) = msh_disk_p4.p1(is_bnd_dg);
msh_winslow = elliptic_smoothing(msh_opt, curvep1, false);
fprintf('\nafter Winslow:\n  nelem = %d\n', size(msh_winslow.t, 2));
plot_linear(msh_winslow, fullfile(figdir, 'mesh_winslow.png'));

fprintf('\nSaved figures to %s/\n', figdir);


function m = boundary_dg_mask(msh)
data = dginit(msh);
[~, edg] = getbndnodes(msh, data);
[ns, d, nt] = size(msh.p1);
m2 = false(ns, nt);
m2(edg) = true;
m = repmat(reshape(m2, ns, 1, nt), 1, d, 1);
end


% ======================================================================
function plot_linear(msh, filename)
figure; clf;
t = double(msh.t) + 1;
triplot(t', msh.p(1,:), msh.p(2,:), 'k-');
axis equal; axis off;
set(gcf, 'Position', [100 100 600 600]);
exportgraphics(gcf, filename, 'Resolution', 200);
end


function loops = debug_extract_loops(msh)
t   = double(msh.t) + 1;
t2t = double(msh.t2t);
nt  = size(t, 2);
edges = zeros(0, 2);
for it = 1:nt
    v1 = t(1, it); v2 = t(2, it); v3 = t(3, it);
    if t2t(1, it) < 0, edges(end+1, :) = [v2, v3]; end %#ok<AGROW>
    if t2t(2, it) < 0, edges(end+1, :) = [v3, v1]; end %#ok<AGROW>
    if t2t(3, it) < 0, edges(end+1, :) = [v1, v2]; end %#ok<AGROW>
end
nedges = size(edges, 1);
used = false(nedges, 1);
loops = {};
for start = 1:nedges
    if used(start), continue; end
    loop = edges(start, 1);
    curr = edges(start, 2);
    used(start) = true;
    steps = 0;
    while curr ~= loop(1) && steps < nedges
        nxt = find(~used & edges(:,1) == curr, 1);
        if isempty(nxt), break; end
        used(nxt) = true;
        loop(end+1) = curr; %#ok<AGROW>
        curr = edges(nxt, 2);
        steps = steps + 1;
    end
    loops{end+1} = loop; %#ok<AGROW>
end
end


function plot_linear_with_cut(msh, cut_path, filename)
figure; clf;
t = double(msh.t) + 1;
triplot(t', msh.p(1,:), msh.p(2,:), 'Color', [0.5 0.5 0.5]);
hold on;
% Draw the cut path (in original indexing, so it lives on side A).
plot(msh.p(1, cut_path), msh.p(2, cut_path), 'r-', 'LineWidth', 2);
plot(msh.p(1, cut_path), msh.p(2, cut_path), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
axis equal; axis off;
set(gcf, 'Position', [100 100 600 600]);
exportgraphics(gcf, filename, 'Resolution', 200);
end
