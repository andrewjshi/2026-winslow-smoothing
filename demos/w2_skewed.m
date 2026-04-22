% w2_skewed.m - Baseline test for the W-modification experiment.
%
% Starts from w2.m's structured 11x11 mesh, applies a smooth nonlinear
% coordinate stretch  xi -> xi + 0.3*xi*(1-xi)*eta  to skew the reference
% geometry (valid, variable Jacobian, not an affine shear). Applies the
% standard w2 boundary perturbation on side 4 and runs unmodified
% Winslow. Goal: verify that reference skew propagates into the smoothed
% physical mesh -- i.e. that there is actually something for the
% W-modification to fix.

porder = 4;
n = 10;
msh = mshsquare(n+1, n+1);
msh = nodealloc(msh, porder);

figdir = 'figures/w2_skewed';
if ~exist(figdir, 'dir'), mkdir(figdir); end

% --- Skew the reference: move one linear vertex to create a single -------
% sliver cluster in the top-right.
moves = [
    0.8, 0.8,  0.72, 0.79;
];
for k = 1:size(moves, 1)
    from = moves(k, 1:2)';
    to   = moves(k, 3:4)';
    [~, ix] = min(vecnorm(double(msh.p) - from, 2, 1));
    msh.p(:, ix) = to;
end
msh = nodealloc(msh, porder);  % rebuild p1 from the modified p

fprintf('Reference mesh skew: min J = %.4g\n', min(geojac(msh), [], 'all'));

% --- Figure 1: skewed reference ------------------------------------------
figure(1); clf;
dgmeshplot_curved(msh, 4, 1, 0);
set(gcf, 'Position', [114 1 560 420]);
exportgraphics(gcf, fullfile(figdir, 'mesh_reference.png'), 'Resolution', 200);

% --- Boundary perturbation on side 4 (same as w2) -----------------------
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

% --- Figure 3: smoothed (baseline) --------------------------------------
figure(3); clf;
dgmeshplot_curved(msh1, 4, 1, 0);
set(gcf, 'Position', [114 1 560 420]);
exportgraphics(gcf, fullfile(figdir, 'mesh_smoothed.png'), 'Resolution', 200);

% --- Quality metrics per element: scaled J + shape distortion eta --------
% eta is the Liu-Joe inverse mean-ratio shape distortion metric,
% eta = |J|_F^2 / (d * |det J|^(2/d)), range [1, inf),
% with eta = 1 for conformal elements (see Shi-Persson 2020, Eq. 2).
[I_ref, eta_ref] = qualmetrics2d(msh);
[I_tan, eta_tan] = qualmetrics2d(msh_tangled);
[I_sm,  eta_sm]  = qualmetrics2d(msh1);

fprintf('\n              reference    tangled     smoothed\n');
fprintf('min  I      %10.4g  %10.4g  %10.4g\n', min(I_ref), min(I_tan), min(I_sm));
fprintf('max  eta    %10.4g  %10.4g  %10.4g\n', max(eta_ref), max(eta_tan), max(eta_sm));
fprintf('mean eta    %10.4g  %10.4g  %10.4g\n', mean(eta_ref), mean(eta_tan), mean(eta_sm));

% Scaled Jacobian histograms
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

% Shape distortion metric histograms (Liu-Joe / Shi-Persson Eq. 2)
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
% Per-element scaled Jacobian I = min J / max J (Fortunato-Persson 2016)
% and Liu-Joe shape distortion eta = |J|_F^2 / (d |det J|^(2/d))
% (Shi-Persson 2020, Eq. 2), both worst-case over Gauss points.
d = size(msh.p, 1);
data = dginit(msh);
p1x = permute(msh.p1(:,1,:), [1,3,2]);
p1y = permute(msh.p1(:,2,:), [1,3,2]);
phiX = permute(data.gfs(:,2,:), [1,3,2]);
phiY = permute(data.gfs(:,3,:), [1,3,2]);
xX = phiX' * p1x;
xY = phiY' * p1x;
yX = phiX' * p1y;
yY = phiY' * p1y;
detJ = xX .* yY - xY .* yX;                                 % (ngauss, nt)
Jfro2 = xX.^2 + xY.^2 + yX.^2 + yY.^2;                      % |J|_F^2
etap = Jfro2 ./ (d * max(abs(detJ), eps).^(2/d));           % |J|_F^2 / (d |detJ|^(2/d))
I   = (min(detJ, [], 1) ./ max(detJ, [], 1))';              % per element
eta = (max(etap, [], 1))';                                   % per element (worst-case)
end
