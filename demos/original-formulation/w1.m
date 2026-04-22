% w1.m - Fortunato-Persson paradigm on the unit disk.
%
% Reference is a linear triangulation of the disk promoted to the true
% curved boundary via mshcurved. The Dirichlet boundary data is the
% straight-sided p1 (boundary nodes sitting on the chords between
% linear vertices), so the initial physical mesh is tangled along the
% curved boundary. Unmodified Winslow untangles it. Saves three
% figures:
%   figures/w1/mesh_reference.png  -- curved reference
%   figures/w1/mesh_tangled.png    -- initial configuration
%   figures/w1/mesh_smoothed.png   -- Winslow result

porder = 4;
msh = mshcircle(1);
msh = nodealloc(msh, porder);
curvep1 = msh.p1;

msh = mshcurved(msh, []);
msh = nodealloc(msh, porder);

figdir = 'figures/w1';
if ~exist(figdir, 'dir'), mkdir(figdir); end

% --- Figure 1: curved reference mesh ------------------------------------
figure(1); clf;
dgmeshplot_curved(msh, 4, 1, 0);
set(gcf, 'Position', [114 1 560 420]);
exportgraphics(gcf, fullfile(figdir, 'mesh_reference.png'), 'Resolution', 200);

% --- Figure 2: tangled initial configuration ----------------------------
msh_tangled = msh;
msh_tangled.p1 = curvep1;
figure(2); clf;
dgmeshplot_curved(msh_tangled, 4, 1, 0);
set(gcf, 'Position', [114 1 560 420]);
exportgraphics(gcf, fullfile(figdir, 'mesh_tangled.png'), 'Resolution', 200);

% --- Run unmodified Winslow ---------------------------------------------
doplot = false;
msh1 = elliptic_smoothing(msh, curvep1, doplot);

% --- Figure 3: final smoothed mesh --------------------------------------
figure(3); clf;
dgmeshplot_curved(msh1, 4, 1, 0);
set(gcf, 'Position', [114 1 560 420]);
exportgraphics(gcf, fullfile(figdir, 'mesh_smoothed.png'), 'Resolution', 200);

fprintf('\nSaved figures to %s/\n', figdir);
