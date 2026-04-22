% w_compare.m - Side-by-side: reference, unmodified Winslow, W-modified.
% Single sliver at top-right. Saves 3 clean mesh figures (no dot overlay).

porder = 4;
n = 10;
msh = mshsquare(n+1, n+1);
msh = nodealloc(msh, porder);
moves = [0.8 0.8 0.72 0.79];
[~, ix] = min(vecnorm(double(msh.p) - moves(1:2)', 2, 1));
msh.p(:, ix) = moves(3:4);
msh = nodealloc(msh, porder);

[~, edg] = getbndnodes(msh, dginit(msh), 4);
xx = msh.p1(:,1,:); yy = msh.p1(:,2,:);
newx = xx; newy = yy;
newx(edg) = xx(edg) + 0.5*sin(pi*yy(edg));
curvep1 = cat(2, newx, newy);

mA = elliptic_smoothing(msh, curvep1, false);
msh.W_h = compute_Wh_nodal(msh, 'equilateral', 'p1');
mB = elliptic_smoothing_W(msh, curvep1, false);

figdir = 'figures/compare_W';
if ~exist(figdir,'dir'), mkdir(figdir); end

figure(1); clf; dgmeshplot_curved(msh, 4, 0, 0);
title('reference: sliver at (0.72, 0.79)');
set(gcf, 'Position', [100 100 600 600]);
exportgraphics(gcf, fullfile(figdir,'1_reference.png'), 'Resolution', 200);

figure(2); clf; dgmeshplot_curved(mA, 4, 0, 0);
title('unmodified Winslow: sliver -> (0.7964, 0.79)');
set(gcf, 'Position', [100 100 600 600]);
exportgraphics(gcf, fullfile(figdir,'2_unmodified.png'), 'Resolution', 200);

figure(3); clf; dgmeshplot_curved(mB, 4, 0, 0);
title('W-modified (P1): sliver -> (0.8029, 0.79)');
set(gcf, 'Position', [100 100 600 600]);
exportgraphics(gcf, fullfile(figdir,'3_W_modified.png'), 'Resolution', 200);

fprintf('\nSaved to %s/\n', figdir);
