% comb_cut_diagnostic.m - Re-build the comb mesh deterministically and
% highlight vertices 17 and 20 (the 2-cut that breaks 3-connectedness).

rng(1);

% Same comb parameters as comb_tangled.m
W       = 0.92;
H_base  = 0.20;
H_top   = 0.70;
margin  = 0.10;
tooth_w = 0.12;
gap_w   = 0.08;
n_teeth = 4;
n_per_edge_density = 10;
n_interior = 60;

% Corner polygon
C = [0, 0; W, 0; W, H_base];
x = W - margin;
for k = 1:n_teeth
    C = [C; x, H_base; x, H_top; x - tooth_w, H_top; x - tooth_w, H_base]; %#ok<AGROW>
    x = x - tooth_w - gap_w;
end
C = [C; 0, H_base];

% Boundary samples
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

% Interior rejection sampling
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

% Plot
figure(1); clf;
triplot(tri, all_pts(:,1), all_pts(:,2), 'Color', [0.6 0.6 0.6]);
hold on;
plot(all_pts(:,1), all_pts(:,2), 'k.', 'MarkerSize', 6);

% Highlight vertices 17 and 20
v17 = all_pts(17, :);
v20 = all_pts(20, :);
plot(v17(1), v17(2), 'ro', 'MarkerSize', 14, 'LineWidth', 2);
plot(v20(1), v20(2), 'ro', 'MarkerSize', 14, 'LineWidth', 2);
text(v17(1) + 0.02, v17(2), '17', 'Color', 'r', 'FontSize', 14, 'FontWeight', 'bold');
text(v20(1) - 0.05, v20(2), '20', 'Color', 'r', 'FontSize', 14, 'FontWeight', 'bold');

% Shade the "pendant" region above the cut (tooth tip)
% tooth 1 is at x in [0.70, 0.82], tip above y=0.60
tip_poly = [0.70 0.60; 0.82 0.60; 0.82 0.70; 0.70 0.70];
fill(tip_poly(:,1), tip_poly(:,2), [1 0.8 0.8], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
uistack(findobj(gca, 'Type', 'line'), 'top');

% Draw the cut line segment between 17 and 20
plot([v17(1) v20(1)], [v17(2) v20(2)], 'r--', 'LineWidth', 2);

axis equal;
xlim([0.60 0.90]); ylim([0.55 0.75]);
title('Zoom on tooth 1: vertices 17, 20 form a 2-cut');
set(gcf, 'Position', [100 100 800 500]);

figdir = 'figures/comb_tangled';
exportgraphics(gcf, fullfile(figdir, 'cut_diagnostic_zoom.png'), 'Resolution', 200);

% Full view
figure(2); clf;
triplot(tri, all_pts(:,1), all_pts(:,2), 'Color', [0.6 0.6 0.6]);
hold on;
plot(all_pts(:,1), all_pts(:,2), 'k.', 'MarkerSize', 4);
plot(v17(1), v17(2), 'ro', 'MarkerSize', 12, 'LineWidth', 2);
plot(v20(1), v20(2), 'ro', 'MarkerSize', 12, 'LineWidth', 2);
text(v17(1) + 0.01, v17(2) + 0.01, '17', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
text(v20(1) - 0.03, v20(2) + 0.01, '20', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
plot([v17(1) v20(1)], [v17(2) v20(2)], 'r--', 'LineWidth', 2);
axis equal; axis off;
set(gcf, 'Position', [100 100 900 600]);
exportgraphics(gcf, fullfile(figdir, 'cut_diagnostic_full.png'), 'Resolution', 200);

fprintf('v17 = (%.3f, %.3f)\n', v17(1), v17(2));
fprintf('v20 = (%.3f, %.3f)\n', v20(1), v20(2));
fprintf('Saved cut_diagnostic_{zoom,full}.png\n');
