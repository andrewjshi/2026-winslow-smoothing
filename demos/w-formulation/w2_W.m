% w2_W.m - w2.m with the W-modification.
%
% Same reference mesh and boundary perturbation as
% demos/original-formulation/w2.m, but the smoother is replaced by
% elliptic_smoothing_W, fed with a nodal W_h built from an ideal
% equilateral triangle via compute_Wh_nodal (default 'vhp' projection,
% area weighting). Saves three figures:
%   figures/w2_W/mesh_reference.png
%   figures/w2_W/mesh_tangled.png
%   figures/w2_W/mesh_smoothed.png

porder = 4;
n = 10;
msh = mshsquare(n+1, n+1);
msh = nodealloc(msh, porder);

% Nodal W_h from equilateral ideal (lumped L2 projection, V_h^p).
msh.W_h = compute_Wh_nodal(msh, 'equilateral');

figdir = 'figures/w2_W';
if ~exist(figdir, 'dir'), mkdir(figdir); end

% --- Figure 1: reference mesh -------------------------------------------
figure(1); clf;
dgmeshplot_curved(msh, 4, 1, 0);
set(gcf, 'Position', [114 1 560 420]);
exportgraphics(gcf, fullfile(figdir, 'mesh_reference.png'), 'Resolution', 200);

% --- Boundary perturbation on side 4 (same as w2.m) ---------------------
[~, edg] = getbndnodes(msh, dginit(msh), 4);
x = msh.p1(:,1,:);
y = msh.p1(:,2,:);
newx = x;
newy = y;
newx(edg) = x(edg) + 0.5*sin(pi*y(edg));
curvep1 = cat(2, newx, newy);

msh_tangled = msh;
msh_tangled.p1 = curvep1;

% --- Figure 2: tangled initial configuration ----------------------------
figure(2); clf;
dgmeshplot_curved(msh_tangled, 4, 1, 0);
set(gcf, 'Position', [114 1 560 420]);
exportgraphics(gcf, fullfile(figdir, 'mesh_tangled.png'), 'Resolution', 200);

% --- Run W-modified Winslow ---------------------------------------------
doplot = false;
msh1 = elliptic_smoothing_W(msh, curvep1, doplot);

% --- Figure 3: final smoothed mesh --------------------------------------
figure(3); clf;
dgmeshplot_curved(msh1, 4, 1, 0);
set(gcf, 'Position', [114 1 560 420]);
exportgraphics(gcf, fullfile(figdir, 'mesh_smoothed.png'), 'Resolution', 200);

fprintf('\nSaved figures to %s/\n', figdir);
