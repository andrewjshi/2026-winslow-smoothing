% comb_2cut_locations.m - Build the finer comb mesh and visualize the
% locations of the four 2-cuts on the constrained Delaunay graph.

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

% --- Build comb boundary -------------------------------------------------
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

% --- Find all 2-cuts using the fast algorithm ---------------------------
I = reshape(tri(:, [1 2 3 1 2 3])', [], 1);
J = reshape(tri(:, [2 3 1 3 1 2])', [], 1);
A = (sparse(I, J, 1, np, np) + sparse(J, I, 1, np, np)) > 0;

G = graph(A);
art_full = articulation_pts(G);
if any(art_full)
    fprintf('Graph has articulation points (1-cuts), not just 2-cuts.\n');
end

% Find all 2-cuts (boundary first, but also interior pairs)
cuts = zeros(0, 2);
for v = 1:n_b
    keep = true(np, 1); keep(v) = false;
    A_sub = A(keep, keep);
    if nnz(A_sub) == 0, continue; end
    G_sub = graph(A_sub);
    if numedges(G_sub) == 0, continue; end
    art_sub = articulation_pts(G_sub);
    sub_to_orig = find(keep);
    for u_sub = find(art_sub)'
        u = sub_to_orig(u_sub);
        cuts(end+1, :) = sort([v, u]); %#ok<AGROW>
    end
end
cuts = unique(cuts, 'rows');

fprintf('\nFound %d 2-cuts:\n', size(cuts, 1));
for c = 1:size(cuts, 1)
    v1 = cuts(c, 1); v2 = cuts(c, 2);
    fprintf('  cut %d: vertices (%d, %d) at (%.4f, %.4f) and (%.4f, %.4f)\n', ...
        c, v1, v2, all_pts(v1,1), all_pts(v1,2), all_pts(v2,1), all_pts(v2,2));
end

% --- Visualize ----------------------------------------------------------
figdir = 'figures/comb';
if ~exist(figdir, 'dir'), mkdir(figdir); end

figure(1); clf;
% Draw mesh
triplot(tri, all_pts(:, 1), all_pts(:, 2), 'Color', [0.6 0.6 0.6]);
hold on;
% Draw boundary points (grey, small)
plot(bnd_pts(:, 1), bnd_pts(:, 2), '.', 'Color', [0.4 0.4 0.4], 'MarkerSize', 4);

% Highlight each 2-cut: draw the two endpoints, plus a connecting line, in red
colors = lines(size(cuts, 1));
for c = 1:size(cuts, 1)
    v1 = cuts(c, 1); v2 = cuts(c, 2);
    p1 = all_pts(v1, :); p2 = all_pts(v2, :);
    plot([p1(1) p2(1)], [p1(2) p2(2)], '-', 'Color', colors(c,:), 'LineWidth', 2);
    plot([p1(1) p2(1)], [p1(2) p2(2)], 'o', ...
        'MarkerFaceColor', colors(c,:), 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 12, 'LineWidth', 1.2);
    % Label
    midp = (p1 + p2) / 2;
    text(midp(1), midp(2) + 0.03, sprintf('cut %d', c), ...
        'Color', colors(c,:), 'FontWeight', 'bold', 'FontSize', 11, ...
        'HorizontalAlignment', 'center');
end

axis equal; axis off;
set(gcf, 'Position', [100 100 1100 600]);
exportgraphics(gcf, fullfile(figdir, '2cuts_locations.png'), 'Resolution', 200);
close;

fprintf('\nSaved figure to %s/2cuts_locations.png\n', figdir);


% =========================================================================
function art = articulation_pts(G)
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
