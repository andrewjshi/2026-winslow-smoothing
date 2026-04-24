% bump.m - Fortunato-Persson 2016 paradigm on the High-Order CFD
% Workshop Gaussian-bump geometry.
%
% Structured rectangular reference mesh of [-1.5, 1.5] x [0, 0.8];
% the bottom boundary (y = 0) is pushed up onto the Gaussian wall
%   y_wall(x) = 0.0625 * exp(-25 x^2)
% The resulting physical high-order mesh is locally tangled in a
% neighbourhood of x = 0 where the bump peak exceeds an element
% layer; Winslow-Picard smoothing untangles. Saves three figures:
%   figures/bump/mesh_reference.png  -- reference (computational) mesh
%   figures/bump/mesh_tangled.png    -- initial configuration, locally tangled
%   figures/bump/mesh_smoothed.png   -- final mesh after Picard iterations

porder = 4;

% Grid resolution: dx = 0.1, dy = 0.04 -> bump peak crosses ~1.5
% element layers at x = 0, giving a visible tangled band around the
% centreline.
m = 31;   % vertices in x
n = 21;   % vertices in y

% Physical rectangle [-1.5, 1.5] x [0, 0.8]
Lx = 3.0;  x0 = -1.5;
Ly = 0.8;

% Bump parameters (HOCFD benchmark)
bump_amp = 0.0625;
bump_k   = 25;

% --- Reference mesh: unit square scaled to the channel rectangle --------
msh = mshsquare(m, n);
msh = nodealloc(msh, porder);
msh.p(1,:)     = Lx * msh.p(1,:)     + x0;
msh.p(2,:)     = Ly * msh.p(2,:);
msh.p1(:,1,:)  = Lx * msh.p1(:,1,:)  + x0;
msh.p1(:,2,:)  = Ly * msh.p1(:,2,:);

figdir = 'figures/bump';
if ~exist(figdir, 'dir'), mkdir(figdir); end

% --- Figure 1: reference mesh -------------------------------------------
figure(1); clf;
dgmeshplot_curved(msh, 4, 0, 0);
set(gcf, 'Position', [114 1 1000 300]);
exportgraphics(gcf, fullfile(figdir, 'mesh_reference.png'), 'Resolution', 200);

% --- Build the bump perturbation on side 1 (y = 0 boundary) -------------
% mshsquare labels: 1 = bottom (y~0), 2 = right, 3 = top, 4 = left.
[~, edg] = getbndnodes(msh, dginit(msh), 1);
x = msh.p1(:,1,:);
y = msh.p1(:,2,:);
newx = x;
newy = y;
newy(edg) = bump_amp * exp(-bump_k * x(edg).^2);
curvep1 = cat(2, newx, newy);

msh_tangled = msh;
msh_tangled.p1 = curvep1;

% --- Figure 2: tangled initial configuration -----------------------------
figure(2); clf;
dgmeshplot_curved(msh_tangled, 4, 0, 0);
set(gcf, 'Position', [114 1 1000 300]);
exportgraphics(gcf, fullfile(figdir, 'mesh_tangled.png'), 'Resolution', 200);

% --- Run Winslow-Picard -------------------------------------------------
doplot = false;
msh1 = elliptic_smoothing(msh, curvep1, doplot);

% --- Figure 3: final smoothed mesh --------------------------------------
figure(3); clf;
dgmeshplot_curved(msh1, 4, 0, 0);
set(gcf, 'Position', [114 1 1000 300]);
exportgraphics(gcf, fullfile(figdir, 'mesh_smoothed.png'), 'Resolution', 200);

fprintf('\nSaved figures to %s/\n', figdir);
