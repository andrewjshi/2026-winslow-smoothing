% w2_skewed.m - same as w2.m, but with a skewed reference mesh.
%
% One linear vertex of the structured 11x11 reference is moved toward a
% neighbour, producing a small sliver cluster in the upper-right of the
% reference. The boundary perturbation and the unmodified Fortunato-
% Persson smoother are otherwise identical to w2.m. Saves three figures:
%   figures/w2_skewed/mesh_reference.png  -- skewed reference
%   figures/w2_skewed/mesh_tangled.png    -- initial configuration
%   figures/w2_skewed/mesh_smoothed.png   -- Winslow result

porder = 4;
n = 10;
msh = mshsquare(n+1, n+1);
msh = nodealloc(msh, porder);

% Move one interior vertex to create a sliver cluster at (0.8, 0.8).
moves = [0.8, 0.8,  0.72, 0.79];
[~, ix] = min(vecnorm(double(msh.p) - moves(1:2)', 2, 1));
msh.p(:, ix) = moves(3:4);
msh = nodealloc(msh, porder);

figdir = 'figures/w2_skewed';
if ~exist(figdir, 'dir'), mkdir(figdir); end

% --- Figure 1: skewed reference mesh ------------------------------------
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

% --- Run unmodified Winslow ---------------------------------------------
doplot = false;
msh1 = elliptic_smoothing(msh, curvep1, doplot);

% --- Figure 3: final smoothed mesh --------------------------------------
figure(3); clf;
dgmeshplot_curved(msh1, 4, 1, 0);
set(gcf, 'Position', [114 1 560 420]);
exportgraphics(gcf, fullfile(figdir, 'mesh_smoothed.png'), 'Resolution', 200);

% --- Quality histograms: scaled Jacobian I and shape distortion eta -----
[I_ref, eta_ref] = qualmetrics2d(msh);
[I_tan, eta_tan] = qualmetrics2d(msh_tangled);
[I_sm,  eta_sm ] = qualmetrics2d(msh1);

fprintf('\n              reference    tangled     smoothed\n');
fprintf('min  I      %10.4g  %10.4g  %10.4g\n', min(I_ref), min(I_tan), min(I_sm));
fprintf('max  eta    %10.4g  %10.4g  %10.4g\n', max(eta_ref), max(eta_tan), max(eta_sm));
fprintf('mean eta    %10.4g  %10.4g  %10.4g\n', mean(eta_ref), mean(eta_tan), mean(eta_sm));

Iall = [I_ref(:); I_tan(:); I_sm(:)];
edges_I = linspace(min(Iall), 1, 40);
configs = {
    'reference', I_ref, 'hist_I_reference.png';
    'tangled',   I_tan, 'hist_I_tangled.png';
    'smoothed',  I_sm,  'hist_I_smoothed.png';
};
for k = 1:3
    figure(3+k); clf;
    histogram(configs{k,2}(:), edges_I, 'FaceColor', [.6, .85, .6]);
    xline(0, 'k--');
    xlabel('scaled Jacobian I'); ylabel('element count');
    title(sprintf('%s: min I = %.4g', configs{k,1}, min(configs{k,2})));
    set(gcf, 'Position', [114 1 560 420]); drawnow;
    exportgraphics(gcf, fullfile(figdir, configs{k,3}), 'Resolution', 200);
end

etaall = [eta_ref(:); eta_tan(:); eta_sm(:)];
edges_eta = linspace(1, quantile(etaall, 0.99), 40);
configs = {
    'reference', eta_ref, 'hist_eta_reference.png';
    'tangled',   eta_tan, 'hist_eta_tangled.png';
    'smoothed',  eta_sm,  'hist_eta_smoothed.png';
};
for k = 1:3
    figure(6+k); clf;
    histogram(configs{k,2}(:), edges_eta, 'FaceColor', [.7, .8, .95]);
    xlabel('shape distortion \eta'); ylabel('element count');
    title(sprintf('%s: mean \\eta = %.4g, max = %.4g', ...
        configs{k,1}, mean(configs{k,2}), max(configs{k,2})));
    set(gcf, 'Position', [114 1 560 420]); drawnow;
    exportgraphics(gcf, fullfile(figdir, configs{k,3}), 'Resolution', 200);
end

fprintf('\nSaved figures to %s/\n', figdir);


function [I, eta] = qualmetrics2d(msh)
d = size(msh.p, 1);
data = dginit(msh);
p1x = permute(msh.p1(:,1,:), [1,3,2]);
p1y = permute(msh.p1(:,2,:), [1,3,2]);
phiX = permute(data.gfs(:,2,:), [1,3,2]);
phiY = permute(data.gfs(:,3,:), [1,3,2]);
xX = phiX' * p1x;  xY = phiY' * p1x;
yX = phiX' * p1y;  yY = phiY' * p1y;
detJ = xX .* yY - xY .* yX;
Jfro2 = xX.^2 + xY.^2 + yX.^2 + yY.^2;
etap = Jfro2 ./ (d * max(abs(detJ), eps).^(2/d));
I   = (min(detJ, [], 1) ./ max(detJ, [], 1))';
eta = (max(etap, [], 1))';
end
