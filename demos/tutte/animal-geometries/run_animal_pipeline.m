function run_animal_pipeline(name, bnd, int_pts, porder, varargin)
% RUN_ANIMAL_PIPELINE  Build a mesh from a closed-polygon boundary and
% interior seed points, diagnose/fix 3-connectedness via interior anchoring,
% and run Tutte + shape-opt against three targets (square, circle, and a
% rectangle matching the domain aspect ratio).
%
% Usage:
%   run_animal_pipeline(name, bnd, int_pts, porder)
%   run_animal_pipeline(..., 'MaxAnchorIter', 10)
%   run_animal_pipeline(..., 'MaxShapeIter', 1000)
%
% Saves figures to figures/<name>/.

opts.MaxAnchorIter = 10;
opts.MaxShapeIter  = 1000;
for k = 1:2:numel(varargin)
    opts.(varargin{k}) = varargin{k+1};
end

figdir = fullfile('figures', name);
if ~exist(figdir, 'dir'), mkdir(figdir); end

n_b = size(bnd, 1);

% --- Interior anchoring loop to restore 3-connectedness ------------------
anchored_int = int_pts;
iter = 0;
while true
    iter = iter + 1;
    msh = build_msh(bnd, anchored_int, porder);
    cuts = find_all_2_cuts(msh, n_b);
    fprintf('[%s] iter %d: %d interior pts, %d 2-cuts\n', ...
        name, iter, size(anchored_int, 1), size(cuts, 1));
    for c = 1:min(size(cuts, 1), 6)
        v1 = cuts(c, 1); v2 = cuts(c, 2);
        fprintf('    cut (%d,%d): (%.3f,%.3f) -- (%.3f,%.3f)\n', ...
            v1, v2, msh.p(1,v1), msh.p(2,v1), msh.p(1,v2), msh.p(2,v2));
    end
    if isempty(cuts) || iter >= opts.MaxAnchorIter, break; end
    % For each 2-cut, find the smaller separated component and place an
    % anchor at its centroid. This is more robust than the midpoint of the
    % cut pair, which for a narrow fin often lands in the main body rather
    % than inside the fin.
    added = 0;
    np_all = size(msh.p, 2);
    t_all  = double(msh.t) + 1;
    I = reshape(t_all([1 2 3 1 2 3], :), [], 1);
    J = reshape(t_all([2 3 1 3 1 2], :), [], 1);
    A_all = (sparse(I, J, 1, np_all, np_all) + sparse(J, I, 1, np_all, np_all)) > 0;
    used_cuts = false(size(cuts, 1), 1);
    for c = 1:size(cuts, 1)
        if used_cuts(c), continue; end
        v1 = cuts(c, 1); v2 = cuts(c, 2);
        keep = true(np_all, 1); keep(v1) = false; keep(v2) = false;
        comps = connected_components(A_all, keep);
        if numel(comps) < 2, continue; end
        [~, smallest] = min(cellfun(@numel, comps));
        verts = comps{smallest};
        candidates = [mean(msh.p(:, verts), 2)'; ...
                      (msh.p(:, v1) + msh.p(:, v2))' / 2];
        domain_center = mean(bnd, 1);
        placed = false;
        for ci = 1:size(candidates, 1)
            anchor = candidates(ci, :);
            % Nudge inward toward domain centroid until >=0.005 from any
            % existing vertex. This handles pointy-corner cuts where the
            % cut-off component is a single boundary vertex, so centroid
            % and midpoint both land ON the boundary.
            for nudge = 0:5
                if ~inpolygon(anchor(1), anchor(2), bnd(:,1), bnd(:,2)), break; end
                existing = [bnd; anchored_int];
                if min(vecnorm(existing - anchor, 2, 2)) >= 0.005
                    anchored_int(end+1, :) = anchor;  %#ok<AGROW>
                    added = added + 1;
                    placed = true;
                    break;
                end
                step = domain_center - anchor;
                step = step / max(norm(step), eps);
                anchor = anchor + 0.015 * step;
            end
            if placed, break; end
        end
        if ~placed, continue; end
        % Mark sibling cuts (sharing one endpoint with (v1,v2) and yielding
        % the same small component) as used to avoid piling up anchors.
        for c2 = c+1:size(cuts, 1)
            u1 = cuts(c2, 1); u2 = cuts(c2, 2);
            if u1 == v1 || u1 == v2 || u2 == v1 || u2 == v2
                used_cuts(c2) = true;
            end
        end
    end
    if added == 0, break; end
end

fprintf('[%s] final: nelem = %d, np = %d, anchors added = %d\n', ...
    name, size(msh.t, 2), size(msh.p, 2), size(anchored_int, 1) - size(int_pts, 1));

% --- Input mesh figure ---------------------------------------------------
plot_mesh(msh, fullfile(figdir, 'mesh_input.png'), bbox_figsize(bnd));

% --- Two targets (rectangle was shown to be much worse on the lake test) -
run_target(msh, msh, porder, figdir, 'square', {'TargetShape', 'square'}, ...
    [500 500], opts.MaxShapeIter, bbox_figsize(bnd), '');
run_target(msh, msh, porder, figdir, 'circle', {'TargetShape', 'circle'}, ...
    [500 500], opts.MaxShapeIter, bbox_figsize(bnd), 'circle');

fprintf('[%s] saved figures to %s/\n\n', name, figdir);
end


% =========================================================================
function msh = build_msh(bnd, int_pts, porder)
all_pts = [bnd; int_pts];
n_b = size(bnd, 1);
ce  = [(1:n_b)', [(2:n_b)'; 1]];
dt  = delaunayTriangulation(all_pts, ce);
in  = isInterior(dt);
tri = dt.ConnectivityList(in, :);
msh = ml2msh(all_pts, tri, {'true'});
msh = mshcurved(msh, []);
msh = nodealloc(msh, porder);
end


function cuts = find_all_2_cuts(msh, n_b)
np = size(msh.p, 2);
t  = double(msh.t) + 1;
I = reshape(t([1 2 3 1 2 3], :), [], 1);
J = reshape(t([2 3 1 3 1 2], :), [], 1);
A = sparse(I, J, 1, np, np);
A = (A + A') > 0;

cuts = zeros(0, 2);
for i = 1:np-1
    if i > n_b, continue; end
    for j = i+1:np
        keep = true(np, 1);
        keep(i) = false; keep(j) = false;
        if ~is_connected_sub(A, keep)
            cuts(end+1, :) = [i, j];  %#ok<AGROW>
        end
    end
end
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
            queue(end+1) = nb;  %#ok<AGROW>
        end
    end
end
c = all(visited(verts));
end


function comps = connected_components(A, keep)
% Returns cell array of component vertex-index lists.
verts = find(keep);
visited = false(size(A, 1), 1);
comps = {};
for s = verts(:)'
    if visited(s), continue; end
    queue = s; visited(s) = true;
    comp = s;
    while ~isempty(queue)
        v = queue(end); queue(end) = [];
        nbrs = find(A(v, :) & keep');
        for nb = nbrs
            if ~visited(nb)
                visited(nb) = true;
                queue(end+1) = nb;  %#ok<AGROW>
                comp(end+1) = nb;   %#ok<AGROW>
            end
        end
    end
    comps{end+1} = comp;  %#ok<AGROW>
end
end


function run_target(msh_in, msh_physical, porder, figdir, name, tutte_args, figsize, max_iter, phys_figsize, slide_mode)
fprintf('  target: %s\n', name);
try
    msh_t = tutte_embedding(msh_in, tutte_args{:});
    msh_t = nodealloc(msh_t, porder);
    plot_mesh(msh_t, fullfile(figdir, sprintf('mesh_tutte_%s.png', name)), figsize);
    shape_args = {'MaxIterations', max_iter};
    if ~isempty(slide_mode)
        shape_args = [shape_args, {'SlideBoundary', slide_mode}];
    end
    [msh_o, info] = optimize_shape(msh_t, shape_args{:});
    msh_o = nodealloc(msh_o, porder);
    fprintf('    F_init = %.3e, F_final = %.3e, iters = %d (slide=%s)\n', ...
        info.F_initial, info.F_final, info.iterations, ...
        ternary(isempty(slide_mode), 'none', slide_mode));
    plot_mesh(msh_o, fullfile(figdir, sprintf('mesh_optimized_%s.png', name)), figsize);

    % Skip Winslow hand-off if shape-opt output is degenerate (else Winslow
    % can hang trying to untangle an unsalvageable input).
    if info.F_final > 1e8
        fprintf('    skipping Winslow (F_final = %.3e too large)\n', info.F_final);
        return;
    end

    curvep1 = msh_o.p1;
    is_bnd_dg = boundary_dg_mask(msh_o);
    curvep1(is_bnd_dg) = msh_physical.p1(is_bnd_dg);
    msh_w = elliptic_smoothing(msh_o, curvep1, false);
    plot_mesh(msh_w, fullfile(figdir, sprintf('mesh_winslow_%s.png', name)), phys_figsize);
catch err
    fprintf('    FAILED: %s\n', err.message);
end
end


function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end


function m = boundary_dg_mask(msh)
data = dginit(msh);
[~, edg] = getbndnodes(msh, data);
[ns, d, nt] = size(msh.p1);
m2 = false(ns, nt);
m2(edg) = true;
m = repmat(reshape(m2, ns, 1, nt), 1, d, 1);
end


function plot_mesh(msh, filename, figsize)
figure; clf;
dgmeshplot_curved(msh, 4, 0, 0);
axis equal; axis off;
set(gcf, 'Position', [100 100 figsize(1) figsize(2)]);
exportgraphics(gcf, filename, 'Resolution', 200);
close;
end


function sz = bbox_figsize(bnd)
xr = max(bnd(:,1)) - min(bnd(:,1));
yr = max(bnd(:,2)) - min(bnd(:,2));
if xr >= yr
    w = 700;
    h = max(250, round(w * yr / xr));
else
    h = 700;
    w = max(250, round(h * xr / yr));
end
sz = [w, h];
end


function sz = rect_figsize(rect_corners)
xr = max(rect_corners(1,:)) - min(rect_corners(1,:));
yr = max(rect_corners(2,:)) - min(rect_corners(2,:));
if xr >= yr
    w = 700;
    h = max(250, round(w * yr / xr));
else
    h = 700;
    w = max(250, round(h * xr / yr));
end
sz = [w, h];
end
