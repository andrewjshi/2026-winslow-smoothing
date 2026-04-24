% w2_tutte.m - Tutte + shape-optimisation pipeline applied to the
% violently-tangled initial configuration from w2.m.
%
% Pipeline:
%   1. Start from w2.m's tangled mesh (mshsquare + boundary perturbation).
%   2. Discard geometry; Tutte on the linear-vertex topology, targeting
%      the unit square.
%   3. Global shape optimisation of interior vertices (min sum_K eta_K^2).
%   4. Promote back to high order via nodealloc.
%
% The w2.m topology is a regular 11x11 grid, so we expect Tutte to
% recover the uniform unit-square grid from scratch -- even though the
% input boundary is curved by sin(pi y) on side 4. Shape-opt should
% have nothing to do on the recovered grid.

porder = 4;
n      = 10;

msh_sq = mshsquare(n+1, n+1);
msh_sq = nodealloc(msh_sq, porder);
[~, edg] = getbndnodes(msh_sq, dginit(msh_sq), 4);
x = msh_sq.p1(:,1,:); y = msh_sq.p1(:,2,:);
newx = x; newy = y;
newx(edg) = x(edg) + 0.5*sin(pi*y(edg));
msh_input = msh_sq;
msh_input.p1 = cat(2, newx, newy);

figdir = 'figures/w2_tutte';
if ~exist(figdir, 'dir'), mkdir(figdir); end

fprintf('input (tangled w2):\n'); print_quality(msh_input);
plot_mesh(msh_input, fullfile(figdir, 'mesh_input.png'));

% --- Tutte (target = unit square, boundary NOT preserved) ---------------
msh_tutte = tutte_embedding(msh_input, 'TargetShape', 'square');
msh_tutte = nodealloc(msh_tutte, porder);

fprintf('\nafter Tutte:\n'); print_quality(msh_tutte);
plot_mesh(msh_tutte, fullfile(figdir, 'mesh_tutte.png'));

% --- Shape optimisation --------------------------------------------------
[msh_opt, info] = optimize_shape(msh_tutte, 'Display', 'iter');
msh_opt = nodealloc(msh_opt, porder);

fprintf('\nafter shape optimisation:\n'); print_quality(msh_opt);
fprintf('  F_init = %.6e, F_final = %.6e, iterations = %d\n', ...
        info.F_initial, info.F_final, info.iterations);

dp = msh_opt.p - msh_tutte.p;
fprintf('  max |dp| = %.4g, mean |dp| = %.4g\n', ...
        max(abs(dp(:))), mean(abs(dp(:))));
plot_mesh(msh_opt, fullfile(figdir, 'mesh_optimized.png'));

fprintf('\nSaved figures to %s/\n', figdir);


function print_quality(msh)
[eta, I] = quality(msh);
fprintf('  nelem = %d\n', numel(eta));
fprintf('  eta   : min = %.4g, mean = %.4g, max = %.4g\n', ...
        min(eta), mean(eta), max(eta));
fprintf('  I     : min = %.4g (I<=0 => inverted)\n', min(I));
end


function [eta, I] = quality(msh)
d = size(msh.p, 1);
data = dginit(msh);
p1x = permute(msh.p1(:,1,:), [1,3,2]);
p1y = permute(msh.p1(:,2,:), [1,3,2]);
phiX = permute(data.gfs(:,2,:), [1,3,2]);
phiY = permute(data.gfs(:,3,:), [1,3,2]);
xX = phiX'*p1x; xY = phiY'*p1x;
yX = phiX'*p1y; yY = phiY'*p1y;
detJ = xX.*yY - xY.*yX;
Jfro2 = xX.^2 + xY.^2 + yX.^2 + yY.^2;
etap = Jfro2 ./ (d * max(abs(detJ), eps).^(2/d));
eta = max(etap, [], 1)';
I = (min(detJ, [], 1) ./ max(detJ, [], 1))';
end


function plot_mesh(msh, filename)
figure; clf;
dgmeshplot_curved(msh, 4, 1, 0);
set(findobj(gca, 'Type', 'line'), 'MarkerSize', 3);
set(gcf, 'Position', [114 1 560 420]);
exportgraphics(gcf, filename, 'Resolution', 200);
end
