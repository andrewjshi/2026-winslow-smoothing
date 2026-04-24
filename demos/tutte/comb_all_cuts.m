% comb_all_cuts.m - Rebuild the comb clean mesh and enumerate every
% 2-vertex cut whose removal disconnects the linear-vertex graph.
% For each cut, report the graph component sizes (so we can spot
% small "pendant" blocks).

rng(1);

W       = 0.92;
H_base  = 0.20;
H_top   = 0.70;
margin  = 0.10;
tooth_w = 0.12;
gap_w   = 0.08;
n_teeth = 4;
n_per_edge_density = 10;
n_interior = 60;

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
bnd_pts = [];
for i = 1:nc
    L = edge_len(i);
    npts = max(1, round(L * n_per_edge_density));
    for k = 0:npts-1
        t = k / npts;
        bnd_pts = [bnd_pts; C(i, :) + t * edge_vec(i, :)]; %#ok<AGROW>
    end
end

int_pts = zeros(0, 2);
while size(int_pts, 1) < n_interior
    pt = [W * rand(), H_top * rand()];
    if inpolygon(pt(1), pt(2), bnd_pts(:,1), bnd_pts(:,2))
        int_pts = [int_pts; pt]; %#ok<AGROW>
    end
end

all_pts = [bnd_pts; int_pts];
n_b = size(bnd_pts, 1);
ce  = [(1:n_b)', [(2:n_b)'; 1]];
dt  = delaunayTriangulation(all_pts, ce);
in  = isInterior(dt);
tri = dt.ConnectivityList(in, :);

np = size(all_pts, 1);
I = reshape(tri(:, [1 2 3 1 2 3])', [], 1);
J = reshape(tri(:, [2 3 1 3 1 2])', [], 1);
A = sparse(I, J, 1, np, np);
A = (A + A') > 0;

fprintf('np = %d, enumerating all 2-cuts...\n', np);

cuts = [];
for i = 1:np-1
    for j = i+1:np
        keep = true(np, 1); keep(i) = false; keep(j) = false;
        comps = connected_components(A, keep);
        if numel(comps) > 1
            pendant_size = min(cellfun(@numel, comps));
            cuts = [cuts; i, j, pendant_size]; %#ok<AGROW>
        end
    end
end

fprintf('\nFound %d 2-cuts:\n', size(cuts, 1));
fprintf('  %-8s %-8s %-20s %-20s %-12s\n', 'v1', 'v2', 'pos v1', 'pos v2', 'pendant_sz');
for k = 1:size(cuts, 1)
    v1 = cuts(k, 1); v2 = cuts(k, 2);
    fprintf('  %-8d %-8d (%6.3f,%6.3f)     (%6.3f,%6.3f)     %d\n', ...
            v1, v2, all_pts(v1, 1), all_pts(v1, 2), ...
            all_pts(v2, 1), all_pts(v2, 2), cuts(k, 3));
end


function comps = connected_components(A, keep)
verts = find(keep);
comps = {};
if isempty(verts), return; end
visited = false(size(A, 1), 1);
while any(keep & ~visited)
    start = find(keep & ~visited, 1);
    stack = start; visited(start) = true;
    this = start;
    while ~isempty(stack)
        v = stack(end); stack(end) = [];
        nbrs = find(A(v, :) & keep' & ~visited');
        for nb = nbrs
            visited(nb) = true;
            stack(end+1) = nb; %#ok<AGROW>
            this(end+1) = nb; %#ok<AGROW>
        end
    end
    comps{end+1} = this; %#ok<AGROW>
end
end
