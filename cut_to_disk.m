function [msh_out, info] = cut_to_disk(msh)
% CUT_TO_DISK  Slit a multiply connected 2D mesh along shortest graph paths
% between each hole and the outer boundary, turning it into a topological
% disk that Tutte's algorithm can handle.
%
% Usage:
%   [msh_out, info] = cut_to_disk(msh)
%
% Input:
%   msh : dgmatlab 2D mesh struct with fields p (d x np), t (3 x nt, 0-based),
%         t2t (3 x nt).  msh.p1 is ignored; downstream code should re-run
%         nodealloc on msh_out.
%
% Output:
%   msh_out : new msh with p, t, t2t (and relabelled boundary indicator)
%             whose linear-vertex domain has one boundary loop. Interior
%             vertices along each cut have been duplicated so triangles on
%             opposite sides of the slit see different copies.
%   info    : struct with fields
%               n_holes  : number of holes the input had
%               cuts     : cell array of cut-path vertex-index sequences
%                          (in the ORIGINAL msh's indexing for the first
%                          cut, subsequent cuts in the progressively
%                          updated indexing)
%
% Algorithm:
%   1. Extract all directed boundary edges (each keeps the mesh on its
%      left) and walk them into oriented loops. Positive signed area
%      identifies the outer loop; holes come out with negative area.
%   2. For each hole, BFS in the linear-vertex adjacency graph from all
%      outer-loop vertices to any hole-loop vertex, producing a shortest
%      edge-path v_0 -> v_1 -> ... -> v_k.
%   3. Mark each triangle face corresponding to an edge of that path as
%      a "blocked" dual-graph edge. BFS on the triangle dual graph splits
%      triangles into two sides A and B of the slit. One side keeps the
%      original path-vertex indices; the other gets duplicates.
%   4. Triangles whose vertices lie on the path but are on side B have
%      those references rewritten to the duplicate vertices. Triangles
%      on side A are unchanged.
%   5. After all holes have been slit, ml2msh is used to rebuild the
%      mesh with correct t2t and a single 'true' boundary label.

info.cuts = {};
loops = extract_boundary_loops_internal(msh);
info.n_holes = numel(loops) - 1;
if info.n_holes <= 0
    msh_out = msh;
    return;
end

msh_out = msh;
for hole_iter = 1:info.n_holes
    loops = extract_boundary_loops_internal(msh_out);
    areas = cellfun(@(L) signed_area_internal(msh_out.p(:, L)), loops);
    [~, outer_idx] = max(areas);
    hole_candidates = setdiff(1:numel(loops), outer_idx);
    hole_idx = hole_candidates(1);

    outer_verts = loops{outer_idx};
    hole_verts  = loops{hole_idx};

    path = shortest_cut_path_internal(msh_out, outer_verts, hole_verts);
    info.cuts{end+1} = path;
    msh_out = slit_along_path_internal(msh_out, path);
end
end


% =======================================================================
function loops = extract_boundary_loops_internal(msh)
% Returns cell array of loops, each a row vector of 1-based vertex
% indices, directed so that the mesh interior is on the LEFT when
% traversing the loop.  Outer loops come out CCW (positive area); holes
% come out CW (negative area).
t   = double(msh.t) + 1;
t2t = double(msh.t2t);
nt  = size(t, 2);

edges = zeros(0, 2);
for it = 1:nt
    v1 = t(1, it); v2 = t(2, it); v3 = t(3, it);
    % Face j in dgmatlab is opposite vertex j. For a CCW triangle, the
    % directed boundary edge that keeps the interior on its left is:
    %   face 1 (opp v1): v2 -> v3
    %   face 2 (opp v2): v3 -> v1
    %   face 3 (opp v3): v1 -> v2
    if t2t(1, it) < 0, edges(end+1, :) = [v2, v3]; end %#ok<AGROW>
    if t2t(2, it) < 0, edges(end+1, :) = [v3, v1]; end %#ok<AGROW>
    if t2t(3, it) < 0, edges(end+1, :) = [v1, v2]; end %#ok<AGROW>
end

nedges = size(edges, 1);
used   = false(nedges, 1);
loops  = {};
for start = 1:nedges
    if used(start), continue; end
    loop = edges(start, 1);
    curr = edges(start, 2);
    used(start) = true;
    while curr ~= loop(1)
        nxt = find(~used & edges(:,1) == curr, 1);
        if isempty(nxt)
            error('extract_boundary_loops: failed to close loop at vertex %d.', curr);
        end
        used(nxt) = true;
        loop(end+1) = curr; %#ok<AGROW>
        curr = edges(nxt, 2);
    end
    loops{end+1} = loop; %#ok<AGROW>
end
end


% =======================================================================
function s = signed_area_internal(p)
% Shoelace signed area of a 2D polygon; CCW -> positive.
n = size(p, 2);
s = 0;
for k = 1:n
    kn = mod(k, n) + 1;
    s = s + p(1, k) * p(2, kn) - p(1, kn) * p(2, k);
end
s = s / 2;
end


% =======================================================================
function path = shortest_cut_path_internal(msh, outer_verts, hole_verts)
% BFS on the linear-vertex adjacency graph. Multi-source from all outer
% vertices; stop at the first hole vertex reached.
np = size(msh.p, 2);
t  = double(msh.t) + 1;
I  = reshape(t([1 2 3 1 2 3], :), [], 1);
J  = reshape(t([2 3 1 3 1 2], :), [], 1);
A  = sparse(I, J, 1, np, np);
A  = (A + A') > 0;

parent = zeros(np, 1);
dist   = inf(np, 1);
dist(outer_verts) = 0;
% Queue stores vertex indices in BFS order.
queue = outer_verts(:)';
for ov = outer_verts
    parent(ov) = -1;  % sentinel: root
end
is_hole = false(np, 1);
is_hole(hole_verts) = true;

target = 0;
while ~isempty(queue)
    v = queue(1); queue(1) = [];
    if is_hole(v)
        target = v;
        break;
    end
    nbrs = find(A(v, :));
    for nb = nbrs
        if isinf(dist(nb))
            dist(nb) = dist(v) + 1;
            parent(nb) = v;
            queue(end+1) = nb; %#ok<AGROW>
        end
    end
end
if target == 0
    error('shortest_cut_path: hole unreachable from outer boundary.');
end

% Trace back from target to root.
path = target;
while parent(path(1)) > 0
    path = [parent(path(1)), path];
end
% Final path goes outer -> hole, including both endpoints.
end


% =======================================================================
function msh_out = slit_along_path_internal(msh, path)
% Duplicate every vertex on `path` and rewrite triangles on the "B side"
% of the path to refer to the duplicates. The A/B partition comes from
% BFS on the triangle dual graph with cut edges treated as blocked.
np_old = size(msh.p, 2);
t   = double(msh.t) + 1;
t2t = double(msh.t2t);
nt  = size(t, 2);
npath = numel(path);

% 1. Identify blocked dual-graph edges (faces corresponding to path edges).
blocked = false(3, nt);
cut_edges = [path(1:end-1)', path(2:end)'];
ncut = size(cut_edges, 1);
% Build a lookup (a,b) -> cut_edge_index, with both orderings.
cut_map = containers.Map('KeyType', 'char', 'ValueType', 'int32');
for k = 1:ncut
    a = cut_edges(k, 1); b = cut_edges(k, 2);
    cut_map(edge_key(a, b)) = int32(k);
    cut_map(edge_key(b, a)) = int32(k);
end
for it = 1:nt
    for j = 1:3
        v = t(:, it);
        v(j) = [];
        if isKey(cut_map, edge_key(v(1), v(2)))
            blocked(j, it) = true;
        end
    end
end

% 2. Geometric classification: each triangle's side is determined by the
%    signed distance from its centroid to the cut polyline. BFS-based
%    labelling fails on annular domains because the dual graph minus
%    cut edges can remain connected (you can walk around the annulus),
%    so we use pure geometry instead.
side = zeros(nt, 1);
for it = 1:nt
    centroid = mean(msh.p(:, t(:, it)), 2);
    side(it) = side_of_polyline(centroid, msh.p(:, path));
end

% 3. Duplicate path vertices.
path_dup_idx = np_old + (1:npath);
new_p = [msh.p, msh.p(:, path)];

% 4. Rewrite triangles on side B to use duplicates where they touch path.
path_lookup = zeros(np_old, 1);
path_lookup(path) = 1:npath;
new_t = t;
for it = 1:nt
    if side(it) ~= 2, continue; end
    for k = 1:3
        idx = t(k, it);
        if idx <= np_old && path_lookup(idx) > 0
            new_t(k, it) = path_dup_idx(path_lookup(idx));
        end
    end
end

% 5. Rebuild the msh via ml2msh, treating all boundary as a single label.
%    Nudge duplicate positions by a tiny amount so ml2msh cannot merge
%    them back into their originals. The shift is numerically negligible
%    and will be overwritten by Tutte anyway.
eps_shift = 1e-9;
new_p(:, path_dup_idx) = new_p(:, path_dup_idx) + eps_shift;
bndexpr = {'true'};
msh_out = ml2msh(new_p', new_t', bndexpr);
msh_out = mshcurved(msh_out, []);
end


% =======================================================================
function k = edge_key(a, b)
k = sprintf('%d_%d', a, b);
end


% =======================================================================
function s = side_of_polyline(q, poly)
% Returns 1 if q is on the LEFT of the oriented polyline `poly` (2 x n),
% returns 2 if on the RIGHT. Uses the nearest-segment sign test.
n = size(poly, 2);
best_d2 = inf;
best_sign = 0;
for k = 1:n-1
    a = poly(:, k);
    b = poly(:, k+1);
    ab = b - a;
    aq = q - a;
    len2 = ab' * ab;
    if len2 == 0, continue; end
    t = (ab' * aq) / len2;
    t = max(0, min(1, t));
    foot = a + t * ab;
    d2 = norm(q - foot)^2;
    if d2 < best_d2
        best_d2 = d2;
        cross_z = ab(1) * aq(2) - ab(2) * aq(1);
        if cross_z > 0
            best_sign = 1;  % left
        else
            best_sign = 2;  % right
        end
    end
end
s = best_sign;
end
