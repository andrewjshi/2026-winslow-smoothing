% comb_3conn_check.m - Build the comb mesh at the new resolution and just
% check 3-connectedness without trying to fix it. Fast: early-exit at the
% first cut found.

rng(1);

W       = 0.92;
H_base  = 0.20;
H_top   = 0.70;
margin  = 0.10;
tooth_w = 0.12;
gap_w   = 0.08;
n_teeth = 4;
n_per_edge_density = 30;
n_interior = 400;

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
fprintf('np = %d (boundary %d, interior %d), nelem = %d\n', ...
    np, n_b, size(int_pts, 1), size(tri, 1));

% Build adjacency
I = reshape(tri(:, [1 2 3 1 2 3])', [], 1);
J = reshape(tri(:, [2 3 1 3 1 2])', [], 1);
A = (sparse(I, J, 1, np, np) + sparse(J, I, 1, np, np)) > 0;

% Check 3-connectedness with early exit
fprintf('searching for first 2-cut (early exit)...\n');
t0 = tic;
found = false;
for i = 1:np-1
    if i > n_b, break; end
    for j = i+1:np
        keep = true(np, 1); keep(i) = false; keep(j) = false;
        if ~is_connected_sub(A, keep)
            fprintf('  2-cut found at vertices %d, %d in %.2fs\n', i, j, toc(t0));
            fprintf('  positions: (%.4f, %.4f) and (%.4f, %.4f)\n', ...
                all_pts(i,1), all_pts(i,2), all_pts(j,1), all_pts(j,2));
            found = true;
            break;
        end
    end
    if found, break; end
end
if ~found
    fprintf('NO 2-cut found (graph is 3-connected) in %.2fs\n', toc(t0));
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
