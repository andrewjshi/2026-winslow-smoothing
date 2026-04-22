% w2.m - Fortunato-Persson 2016 paradigm, replication and figure generation.
%
% Structured reference mesh, violent boundary perturbation on one side,
% Winslow-Picard smoothing to untangle. Saves three figures:
%   figures/mesh_reference.png  -- reference (computational) mesh
%   figures/mesh_tangled.png    -- initial configuration, heavy tangling
%   figures/mesh_smoothed.png   -- final mesh after Picard iterations

porder = 4;
n = 10;
msh = mshsquare(n+1, n+1);
msh = nodealloc(msh, porder);

% Output directory for figures
figdir = 'figures/w2';
if ~exist(figdir, 'dir'), mkdir(figdir); end

% --- Figure 1: reference mesh -------------------------------------------
figure(1); clf;
dgmeshplot_curved(msh, 4, 1, 0);
set(gcf, 'Position', [114 1 560 420]);
exportgraphics(gcf, fullfile(figdir, 'mesh_reference.png'), 'Resolution', 200);

% --- Build the perturbation: displace nodes on boundary 4 -----------------
% Grid spacing is Delta = 0.1; amplitude 0.5 cuts through roughly five
% element layers, producing heavy tangling near the perturbed boundary.
[~, edg] = getbndnodes(msh, dginit(msh), 4);
x = msh.p1(:,1,:);
y = msh.p1(:,2,:);
newx = x;
newy = y;
newx(edg) = x(edg) + 0.5*sin(pi*y(edg));
curvep1 = cat(2, newx, newy);

msh_tangled = msh;
msh_tangled.p1 = curvep1;

% --- Figure 2: tangled initial configuration -----------------------------
figure(2); clf;
dgmeshplot_curved(msh_tangled, 4, 1, 0);
set(gcf, 'Position', [114 1 560 420]);
exportgraphics(gcf, fullfile(figdir, 'mesh_tangled.png'), 'Resolution', 200);

% --- Run Winslow-Picard ---------------------------------------------------
doplot = false;  % suppress intermediate iteration plots
msh1 = elliptic_smoothing(msh, curvep1, doplot);

% --- Figure 3: final smoothed mesh ---------------------------------------
figure(3); clf;
dgmeshplot_curved(msh1, 4, 1, 0);
set(gcf, 'Position', [114 1 560 420]);
exportgraphics(gcf, fullfile(figdir, 'mesh_smoothed.png'), 'Resolution', 200);

fprintf('\nSaved figures to %s/\n', figdir);
