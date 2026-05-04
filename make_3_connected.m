function [int_pts_out, info] = make_3_connected(bnd, int_pts, varargin)
% MAKE_3_CONNECTED  Add interior anchor points to a 2D simply-connected
% polygonal domain so that the constrained Delaunay triangulation's
% linear-vertex graph becomes 3-connected. This is the precondition for
% Tutte's theorem to produce a valid embedding.
%
% Inputs:
%   bnd      N_b x 2 array, boundary polygon vertices in CCW order.
%   int_pts  N_i x 2 array, initial interior points (may be empty).
%
% Optional name/value pairs:
%   'MaxIter'        Max anchoring iterations (default 20).
%   'Verbose'        Print per-iter diagnostics (default true).
%   'MinSpacing'     Minimum spacing between any new anchor and existing
%                    points (default 0.005). Anchors closer than this are
%                    nudged inward; if no nudge succeeds the cut is skipped
%                    this iteration.
%   'NudgeStep'      Step size when nudging an anchor toward the domain
%                    centroid (default 0.015). Used when an anchor candidate
%                    lands on the boundary or too close to existing points.
%   'BoundaryFirst'  If true (default), only consider 2-cuts where at least
%                    one vertex is on the boundary. This is a 2x speedup and
%                    the typical case --- 2-cuts in CDTs of polygonal
%                    domains almost always include a boundary vertex.
%
% Outputs:
%   int_pts_out   Updated interior points (>= input size).
%   info          Struct: .iters, .added, .final_cuts.

opts.MaxIter       = 20;
opts.Verbose       = true;
opts.MinSpacing    = 0.005;
opts.NudgeStep     = 0.015;
opts.BoundaryFirst = true;
for k = 1:2:numel(varargin)
    opts.(varargin{k}) = varargin{k+1};
end

n_b           = size(bnd, 1);
int_pts_out   = int_pts;
domain_center = mean(bnd, 1);

iter = 0;
added = 0;
last_cuts = [];

while true
    iter = iter + 1;

    % Constrained Delaunay on (boundary, current interior points).
    all_pts = [bnd; int_pts_out];
    np = size(all_pts, 1);
    ce = [(1:n_b)', [(2:n_b)'; 1]];
    dt = delaunayTriangulation(all_pts, ce);
    in = isInterior(dt);
    tri = dt.ConnectivityList(in, :);

    % Linear-vertex graph adjacency.
    I = reshape(tri(:, [1 2 3 1 2 3])', [], 1);
    J = reshape(tri(:, [2 3 1 3 1 2])', [], 1);
    A = (sparse(I, J, 1, np, np) + sparse(J, I, 1, np, np)) > 0;

    % Find all 2-cuts via per-vertex biconnected-components decomposition:
    % for each candidate v, the articulation points of G - v paired with v
    % are exactly the 2-cuts containing v. biconncomp is O(V+E), so total
    % cost is O(V * (V+E)) instead of the O(V^2 * (V+E)) of brute force.
    cuts = find_all_2_cuts(A, n_b, opts.BoundaryFirst);
    last_cuts = cuts;

    if opts.Verbose
        fprintf('  make_3_connected iter %d: %d interior pts, %d 2-cuts\n', ...
            iter, size(int_pts_out, 1), size(cuts, 1));
    end

    if isempty(cuts), break; end
    if iter >= opts.MaxIter
        if opts.Verbose
            fprintf('  make_3_connected: hit MaxIter=%d, %d cuts remain\n', ...
                opts.MaxIter, size(cuts, 1));
        end
        break;
    end

    % Add an anchor for each cut (skipping siblings of cuts already handled).
    used = false(size(cuts, 1), 1);
    iter_added = 0;
    for c = 1:size(cuts, 1)
        if used(c), continue; end
        v1 = cuts(c, 1); v2 = cuts(c, 2);
        keep = true(np, 1); keep(v1) = false; keep(v2) = false;
        comps = connected_components(A, keep);
        if numel(comps) < 2, continue; end
        [~, smallest] = min(cellfun(@numel, comps));
        verts = comps{smallest};
        candidates = [
            mean(all_pts(verts, :), 1);
            (all_pts(v1, :) + all_pts(v2, :)) / 2
        ];
        placed = false;
        for ci = 1:size(candidates, 1)
            anchor = candidates(ci, :);
            for nudge = 0:5
                if ~inpolygon(anchor(1), anchor(2), bnd(:,1), bnd(:,2))
                    break;
                end
                if min(vecnorm(all_pts - anchor, 2, 2)) >= opts.MinSpacing
                    int_pts_out = [int_pts_out; anchor]; %#ok<AGROW>
                    iter_added = iter_added + 1;
                    placed = true;
                    break;
                end
                step = domain_center - anchor;
                step = step / max(norm(step), eps);
                anchor = anchor + opts.NudgeStep * step;
            end
            if placed, break; end
        end
        if ~placed, continue; end
        % Mark sibling cuts (sharing one endpoint) as handled.
        for c2 = c+1:size(cuts, 1)
            u1 = cuts(c2, 1); u2 = cuts(c2, 2);
            if any([u1, u2] == v1) || any([u1, u2] == v2)
                used(c2) = true;
            end
        end
    end
    added = added + iter_added;

    if iter_added == 0
        if opts.Verbose
            fprintf('  make_3_connected: stalled (no anchors placed); %d cuts remain\n', ...
                size(cuts, 1));
        end
        break;
    end
end

info.iters      = iter;
info.added      = added;
info.final_cuts = size(last_cuts, 1);
end


% =========================================================================
function cuts = find_all_2_cuts(A, n_b, boundary_first)
% Find all pairs of vertices (i, j) such that removing both from the
% linear-vertex graph disconnects it. Uses MATLAB's biconncomp:
%
%   - If the full graph has articulation points (1-cuts), those are
%     symptoms of a graph that's not even 2-connected. Emit a 2-cut for
%     each: (articulation_pt, one vertex from the smallest disconnected
%     component) so the anchor-placement loop can repair it.
%
%   - If the graph IS 2-connected, then a 2-cut (v, u) corresponds to: v
%     plus an articulation point u of G - v. We iterate v over candidate
%     vertices, run biconncomp on G - v, and read off articulation points.
%
% Total cost: O(n * (V + E)) instead of O(n^2 * (V + E)).
n = size(A, 1);
G = graph(A);

% Phase 1: articulation points (1-cuts) of the full graph.
art_full = articulation_points(G);
if any(art_full)
    cuts = zeros(0, 2);
    for v = find(art_full)'
        keep = true(n, 1); keep(v) = false;
        comps = connected_components(A, keep);
        if numel(comps) < 2, continue; end
        [~, ix] = min(cellfun(@numel, comps));
        u = comps{ix}(1);
        cuts(end+1, :) = sort([v, u]); %#ok<AGROW>
    end
    return;
end

% Phase 2: graph is biconnected; find true 2-cuts.
v_range = 1:n;
if boundary_first
    v_range = 1:n_b;
end

cuts = zeros(0, 2);
for v = v_range
    keep = true(n, 1); keep(v) = false;
    A_sub = A(keep, keep);
    if nnz(A_sub) == 0, continue; end
    G_sub = graph(A_sub);
    if numedges(G_sub) == 0, continue; end
    art_sub = articulation_points(G_sub);
    sub_to_orig = find(keep);
    for u_sub = find(art_sub)'
        u = sub_to_orig(u_sub);
        cuts(end+1, :) = sort([v, u]); %#ok<AGROW>
    end
end
if ~isempty(cuts)
    cuts = unique(cuts, 'rows');
end
end


function art = articulation_points(G)
% Articulation points = vertices appearing in 2 or more biconnected
% components. biconncomp returns a per-edge bcc index; a vertex is an
% articulation point iff its incident edges span 2+ bcc indices.
n = numnodes(G);
art = false(n, 1);
if numedges(G) == 0, return; end
bins = biconncomp(G);
[src, tgt] = findedge(G);
for v = 1:n
    my = (src == v) | (tgt == v);
    if ~any(my), continue; end
    if numel(unique(bins(my))) >= 2
        art(v) = true;
    end
end
end


function comps = connected_components(A, keep)
% Cell array of connected component vertex-index lists in the subgraph
% induced by keep.
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
                queue(end+1) = nb; %#ok<AGROW>
                comp(end+1) = nb;  %#ok<AGROW>
            end
        end
    end
    comps{end+1} = comp; %#ok<AGROW>
end
end
