% comb_preview.m - Draw the proposed comb domain and a straight-sided
% Delaunay mesh on it. No pipeline run yet; this is for shape approval.

% --- Comb parameters -----------------------------------------------------
W           = 0.92;   % total width
H_base      = 0.20;   % height of the comb's back/base
H_top       = 0.70;   % top of each tooth (tooth height = 0.50)
margin      = 0.10;   % horizontal margin at each end
tooth_w     = 0.12;
gap_w       = 0.08;
n_teeth     = 4;
n_per_side  = 3;      % boundary points per side-length unit of 0.1
n_interior  = 60;

% --- Build corner list (CCW) --------------------------------------------
C = [0, 0; W, 0; W, H_base];
x = W - margin;
for k = 1:n_teeth
    C = [C; x, H_base];              %#ok<AGROW>  inner base point
    C = [C; x, H_top];               %#ok<AGROW>  right side top
    C = [C; x - tooth_w, H_top];     %#ok<AGROW>  left side top
    C = [C; x - tooth_w, H_base];    %#ok<AGROW>  left side base
    x = x - tooth_w - gap_w;
end
C = [C; 0, H_base];     % left margin top
% (last edge closes back to (0, 0))
fprintf('Comb has %d corners.\n', size(C, 1));

% --- Sample boundary points along the polygon edges -----------------------
nc = size(C, 1);
edge_vec = C([2:nc 1], :) - C;
edge_len = vecnorm(edge_vec, 2, 2);
bnd_pts = [];
for i = 1:nc
    L = edge_len(i);
    npts = max(1, round(L * 10 * n_per_side));
    for k = 0:npts-1
        t = k / npts;
        bnd_pts = [bnd_pts; C(i, :) + t * edge_vec(i, :)]; %#ok<AGROW>
    end
end
fprintf('Boundary samples: %d\n', size(bnd_pts, 1));

% --- Interior: rejection into polygon ---------------------------------------
int_pts = zeros(0, 2);
while size(int_pts, 1) < n_interior
    pt = [W * rand(), H_top * rand()];
    if inpolygon(pt(1), pt(2), bnd_pts(:,1), bnd_pts(:,2))
        int_pts = [int_pts; pt]; %#ok<AGROW>
    end
end

all_pts = [bnd_pts; int_pts];

% --- Constrained Delaunay --------------------------------------------------
n_b = size(bnd_pts, 1);
ce  = [(1:n_b)', [(2:n_b)'; 1]];
dt  = delaunayTriangulation(all_pts, ce);
in  = isInterior(dt);
tri = dt.ConnectivityList(in, :);

fprintf('Elements: %d\n', size(tri, 1));

% --- Plot ----------------------------------------------------------------
figdir = 'figures/comb_tangled';
if ~exist(figdir, 'dir'), mkdir(figdir); end

figure(1); clf;
triplot(tri, all_pts(:, 1), all_pts(:, 2), 'k-');
hold on;
plot(bnd_pts(:, 1), bnd_pts(:, 2), 'r.', 'MarkerSize', 8);
axis equal; axis off;
set(gcf, 'Position', [100 100 800 400]);
exportgraphics(gcf, fullfile(figdir, 'comb_mesh.png'), 'Resolution', 200);

figure(2); clf;
fill(C([1:end 1], 1), C([1:end 1], 2), [.85 .95 .85]);
axis equal; axis off;
set(gcf, 'Position', [100 100 800 400]);
exportgraphics(gcf, fullfile(figdir, 'comb_boundary.png'), 'Resolution', 200);

fprintf('Saved figures to %s/\n', figdir);
