% comb.m - Comb-shaped input through the Tutte pipeline, tested against
% both unit-square and unit-circle computational domains.
%
% No interior tangling --- the comb's boundary alone (16 sharp corners,
% four narrow teeth) is a non-trivial untangling problem on its own.
% Mesh resolution is set so each tooth has roughly 3 layers of elements
% across its width.

rng(1);
porder = 4;

% --- Comb parameters -----------------------------------------------------
% Very mild comb: 3 teeth, aspect ~ 1:1 (square nubs of width 0.18). The
% aggressive 4-teeth 4:1 version inverts at beta=0.3; the mild 3-teeth 2:1
% version inverts at beta=0.5; this one tests whether the pipeline can
% complete to beta=1 at all on comb-like geometries.
W       = 0.92;
H_base  = 0.20;
H_top   = 0.38;     % tooth height = H_top - H_base = 0.18 (~1:1 aspect)
margin  = 0.10;
tooth_w = 0.18;
gap_w   = 0.10;
n_teeth = 3;
n_per_edge_density = 25;
n_interior = 350;

% --- Corners of the comb (CCW) ------------------------------------------
C = [0, 0; W, 0; W, H_base];
x = W - margin;
for k = 1:n_teeth
    C = [C; x, H_base; x, H_top; x - tooth_w, H_top; x - tooth_w, H_base]; %#ok<AGROW>
    x = x - tooth_w - gap_w;
end
C = [C; 0, H_base];

% --- Boundary sample points ---------------------------------------------
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

% --- Interior points: rejection sampling --------------------------------
int_pts = zeros(0, 2);
while size(int_pts, 1) < n_interior
    pt = [W * rand(), H_top * rand()];
    if inpolygon(pt(1), pt(2), bnd_pts(:,1), bnd_pts(:,2))
        int_pts = [int_pts; pt]; %#ok<AGROW>
    end
end

% --- Repair 3-connectedness of the linear-vertex graph ------------------
% Tutte's theorem requires 3-connectedness; comb-style domains with thin
% teeth can have constrained-Delaunay 2-cuts at the tooth bases. The helper
% adds interior anchor points until no 2-cut remains.
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

figdir = 'figures/comb';
if ~exist(figdir, 'dir'), mkdir(figdir); end

fprintf('clean comb mesh:\n'); print_quality(msh_clean);
% (3-connectedness already established by make_3_connected above; the
% O(V^2) brute-force diagnostic that lived here was prohibitive at this
% mesh resolution.)
plot_mesh(msh_clean, fullfile(figdir, 'mesh_clean.png'));

% --- Pipeline (circle target only; square target is degenerate because
%     the comb's 16 sharp boundary corners don't fit onto the square's 4
%     corners --- Tutte produces a configuration with eta ~ 1e6).
for target_shape = {'circle'}
    ts = target_shape{1};
    fprintf('\n=== target: %s ===\n', ts);

    msh_t = tutte_embedding(msh_clean, 'TargetShape', ts);
    msh_t = nodealloc(msh_t, porder);
    fprintf('after Tutte:\n'); print_quality(msh_t);
    plot_mesh(msh_t, fullfile(figdir, sprintf('mesh_tutte_%s.png', ts)));

    [msh_o, info] = optimize_shape(msh_t);
    msh_o = nodealloc(msh_o, porder);
    fprintf('after shape-opt:\n'); print_quality(msh_o);
    fprintf('  F_init = %.4e, F_final = %.4e, iterations = %d\n', ...
            info.F_initial, info.F_final, info.iterations);
    plot_mesh(msh_o, fullfile(figdir, sprintf('mesh_optimized_%s.png', ts)));

    % Winslow with beta-homotopy continuation (Fortunato thesis eq. 2.26):
    % step the boundary gradually from the disk reference to the comb
    % target in n_steps stages, each a small, well-conditioned sub-problem.
    try
        n_steps = 10;
        ws_opts = struct('alpha', 0.5, 'maxiter', 100, 'tol', 1e-5);
        [msh_w, hinfo] = winslow_homotopy(msh_o, msh_clean, n_steps, ws_opts);
        fprintf('homotopy: %d/%d stages converged. completed=%d.\n', ...
                hinfo.converged_steps, n_steps, hinfo.completed);
        fprintf('after Winslow:\n'); print_quality(msh_w);
        plot_mesh(msh_w, fullfile(figdir, sprintf('mesh_winslow_%s.png', ts)));
    catch ME
        fprintf('Winslow FAILED for target=%s: %s\n', ts, ME.message);
    end
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
set(gcf, 'Position', [114 1 800 500]);
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
queue = verts(1); visited(queue) = true;
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
