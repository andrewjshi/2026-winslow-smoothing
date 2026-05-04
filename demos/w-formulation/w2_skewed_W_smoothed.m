% w2_skewed_W_smoothed.m -- Sweep of Laplacian-smoothed W_h on the
% skewed-reference test. Same setup as w2_skewed_W.m: skewed reference,
% boundary perturbation on side 4. We replace the unsmoothed P1
% projection of W_h with compute_Wh_smoothed(msh, 'equilateral', n_iter)
% for several n_iter values, run the W-modified Picard solver, and
% report whether the iteration converges and what the final mesh
% quality looks like.
%
% Hypothesis: there is a sweet spot where smoothing is enough to make
% Picard converge but not so much that W_h collapses toward constant
% (which would make the W-modification revert to plain Winslow).

porder = 4;
n = 10;
sweep = [0, 1, 2, 5, 10, 20, 50];

% --- Build the skewed reference mesh once (shared by all variants) ------
msh = mshsquare(n+1, n+1);
msh = nodealloc(msh, porder);
moves = [0.8, 0.8,  0.72, 0.79];
[~, ix] = min(vecnorm(double(msh.p) - moves(1:2)', 2, 1));
msh.p(:, ix) = moves(3:4);
msh = nodealloc(msh, porder);

% --- Boundary perturbation ----------------------------------------------
[~, edg] = getbndnodes(msh, dginit(msh), 4);
x = msh.p1(:,1,:); y = msh.p1(:,2,:);
newx = x; newy = y;
newx(edg) = x(edg) + 0.5*sin(pi*y(edg));
curvep1 = cat(2, newx, newy);

msh_tangled = msh;
msh_tangled.p1 = curvep1;

[I_ref, eta_ref] = qualmetrics2d(msh);
[I_tan, eta_tan] = qualmetrics2d(msh_tangled);
fprintf('reference:  min I = %.4g, max eta = %.4g\n', min(I_ref), max(eta_ref));
fprintf('tangled  :  min I = %.4g, max eta = %.4g\n', min(I_tan), max(eta_tan));

figdir_root = 'figures/w2_skewed_W_smoothed';
if ~exist(figdir_root, 'dir'), mkdir(figdir_root); end

% --- Sweep over smoothing iterations ------------------------------------
results = struct('n_iter', {}, 'converged', {}, 'min_I', {}, 'max_eta', {}, ...
                 'mean_eta', {}, 'gradW_norm', {});

for k = 1:numel(sweep)
    n_iter = sweep(k);
    fprintf('\n========== n_iter = %d ==========\n', n_iter);

    msh_k = msh;
    msh_k.W_h = compute_Wh_smoothed(msh_k, 'equilateral', n_iter);

    % Crude grad(W_h) magnitude: per-element finite difference between
    % the W values at the linear vertices, summed component-wise.
    gradW = estimate_gradW(msh_k);

    converged = true;
    try
        msh_out = elliptic_smoothing_W(msh_k, curvep1, false);
    catch ME
        fprintf('  elliptic_smoothing_W FAILED: %s\n', ME.message);
        msh_out = msh_tangled;
        converged = false;
    end

    [I_sm, eta_sm] = qualmetrics2d(msh_out);
    has_inversions = min(I_sm) <= 0;

    fprintf('  smoothed: min I = %.4g, max eta = %.4g, mean eta = %.4g\n', ...
            min(I_sm), max(eta_sm), mean(eta_sm));
    fprintf('  |grad W_h| (sum_K |W_v_max - W_v_min|) = %.4g\n', gradW);
    if has_inversions
        fprintf('  *** INVERTED ELEMENTS PRESENT ***\n');
    end

    results(k).n_iter     = n_iter;
    results(k).converged  = converged && ~has_inversions;
    results(k).min_I      = min(I_sm);
    results(k).max_eta    = max(eta_sm);
    results(k).mean_eta   = mean(eta_sm);
    results(k).gradW_norm = gradW;

    figdir = fullfile(figdir_root, sprintf('n_iter_%d', n_iter));
    if ~exist(figdir, 'dir'), mkdir(figdir); end
    figure(100); clf;
    dgmeshplot_curved(msh_out, 4, 1, 0);
    set(gcf, 'Position', [114 1 560 420]);
    title(sprintf('n\\_iter = %d, min I = %.3g, max \\eta = %.3g', ...
        n_iter, min(I_sm), max(eta_sm)));
    exportgraphics(gcf, fullfile(figdir, 'mesh_smoothed.png'), 'Resolution', 200);
end

% --- Summary table ------------------------------------------------------
fprintf('\n========== summary ==========\n');
fprintf('%6s  %8s  %10s  %10s  %10s  %12s\n', ...
        'n_iter','clean?','min I','max eta','mean eta','|grad W_h|');
for k = 1:numel(results)
    r = results(k);
    fprintf('%6d  %8s  %10.4g  %10.4g  %10.4g  %12.4g\n', ...
            r.n_iter, ternary(r.converged, 'YES', 'no'), ...
            r.min_I, r.max_eta, r.mean_eta, r.gradW_norm);
end

fprintf('\nSaved figures under %s/\n', figdir_root);


% =========================================================================
function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function gradW = estimate_gradW(msh)
% Sum over elements of |max(W_v) - min(W_v)| over the element's corners,
% component-wise. A coarse proxy for the magnitude of grad(W_h).
nt = size(msh.t, 2);
ns = size(msh.s, 1);
W = msh.W_h;          % (ns, d, d, nt)
gradW = 0;
for it = 1:nt
    Wel = squeeze(W(:, :, :, it));   % (ns, d, d)
    for a = 1:size(Wel, 2)
        for b = 1:size(Wel, 3)
            v = Wel(:, a, b);
            gradW = gradW + (max(v) - min(v));
        end
    end
end
end

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
