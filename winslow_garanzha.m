function [msh_out, info] = winslow_garanzha(msh_ref, bnd_target, varargin)
% WINSLOW_GARANZHA  Foldover-free linear-vertex mesh smoothing via the
% Garanzha 2021 penalty method (Escobar-regularized barrier).
%
% Treats msh_ref as the reference (computational) mesh and computes a
% piecewise-affine map from it to a physical mesh with linear-vertex
% boundary fixed at bnd_target. Minimizes
%
%   E(u) = sum_t [f_eps(J_t) + lambda * g_eps(J_t)] * vol(T_t_ref)
%
% where J_t is the Jacobian of the affine map ref->phys for triangle t,
% with the regularized denominator
%
%   chi(D, eps) = (D + sqrt(eps^2 + D^2)) / 2
%   f_eps(J)    = tr(J^T J) / chi(det J, eps)         (Frobenius / det)
%   g_eps(J)    = (det^2 J + 1) / chi(det J, eps)     (area term)
%
% Outer iteration drives eps -> 0; for det J > 0 this recovers Winslow's
% original variational form, but the regularization keeps the energy and
% its gradient finite for det J <= 0 so inverted initial states untangle.
%
% Inputs:
%   msh_ref     reference mesh (provides connectivity msh.t and reference
%               linear-vertex positions msh.p; msh.t2t for boundary mask)
%   bnd_target  2 x N_b array, target physical positions for boundary
%               linear vertices (in the order returned by find(boundary_mask))
%
% Optional name/value pairs:
%   'Lambda'      area-vs-shape weight (default 1.0)
%   'MaxOuter'    max epsilon outer iterations (default 30)
%   'MaxInner'    max BFGS iterations per outer (default 200)
%   'GradTol'     gradient norm tolerance for inner BFGS (default 1e-7)
%   'EpsInit'     initial epsilon (default 1.0)
%   'EpsTarget'   stop when eps drops below this AND mesh valid (default 1e-9)
%   'Display'     'none' | 'iter' (default 'iter')
%
% Output:
%   msh_out   copy of msh_ref with linear-vertex positions msh.p replaced
%             by the optimized physical positions. High-order DG nodes
%             (msh.p1) are NOT touched; caller should run nodealloc.
%   info      struct: .min_detJ_per_outer, .energy_per_outer, .eps_history,
%                     .converged

opts.Lambda    = 1.0;
opts.MaxOuter  = 30;
opts.MaxInner  = 200;
opts.GradTol   = 1e-7;
opts.EpsInit   = 1.0;
opts.EpsTarget = 1e-9;
opts.Display   = 'iter';
for k = 1:2:numel(varargin)
    opts.(varargin{k}) = varargin{k+1};
end
verbose = strcmpi(opts.Display, 'iter');

p_ref = double(msh_ref.p);
t     = double(msh_ref.t) + 1;
nt    = size(t, 2);
np    = size(p_ref, 2);

is_bnd = boundary_mask(msh_ref);
bnd_v  = find(is_bnd);
int_v  = find(~is_bnd);
nint   = numel(int_v);
nbnd   = numel(bnd_v);

if size(bnd_target, 2) ~= nbnd
    error('winslow_garanzha: bnd_target has %d cols, expected %d.', ...
        size(bnd_target, 2), nbnd);
end

% Initial physical = reference, with boundary replaced by target
p_phys = p_ref;
p_phys(:, bnd_v) = bnd_target;

% Precompute reference edge matrices and Z = (e3xe3 lift) * Sinv per
% triangle, plus reference volumes. For each triangle:
%   S_t (2x2)  = [x1-x0, x2-x0] from REFERENCE positions
%   Z_t (3x2)  = [-1 -1; 1 0; 0 1] * inv(S_t)
%   J_t (2x2)  = U_t (2x3) * Z_t  with U_t = [u0 u1 u2] (current physical)
% so the per-triangle Jacobian is linear in physical vertex positions.
Z       = zeros(3, 2, nt);
volref  = zeros(nt, 1);
for it = 1:nt
    v = t(:, it);
    e1 = p_ref(:, v(2)) - p_ref(:, v(1));
    e2 = p_ref(:, v(3)) - p_ref(:, v(1));
    detS = e1(1)*e2(2) - e1(2)*e2(1);
    if detS <= 0
        warning('winslow_garanzha: reference triangle %d has det <= 0 (%g)', it, detS);
    end
    Sinv = [ e2(2), -e2(1); -e1(2), e1(1) ] / detS;
    Z(:, :, it) = [-1, -1; 1, 0; 0, 1] * Sinv;
    volref(it) = abs(detS) / 2;
end

% Outer epsilon iteration
info.min_detJ_per_outer = nan(opts.MaxOuter, 1);
info.energy_per_outer   = nan(opts.MaxOuter, 1);
info.eps_history        = nan(opts.MaxOuter, 1);
info.converged          = false;

eps_k = opts.EpsInit;

for outer = 1:opts.MaxOuter
    info.eps_history(outer) = eps_k;
    if verbose
        fprintf('  garanzha outer %2d  eps=%.3e ', outer, eps_k);
    end

    % Inner BFGS at current eps over interior positions
    fun = @(x) energy_and_grad(x, p_phys, int_v, t, Z, volref, opts.Lambda, eps_k);
    x0  = reshape(p_phys(:, int_v), [], 1);
    [x_opt, F_opt] = bfgs_minimize(fun, x0, opts.MaxInner, opts.GradTol, false);
    p_phys(:, int_v) = reshape(x_opt, 2, nint);

    % Diagnostics
    [minDJ, ~] = min_det_J(p_phys, t, Z);
    info.min_detJ_per_outer(outer) = minDJ;
    info.energy_per_outer(outer)   = F_opt;

    if verbose
        fprintf('F=%.4e  minDJ=%+.3e\n', F_opt, minDJ);
    end

    % Convergence: valid mesh AND eps below target
    if minDJ > 0 && eps_k <= opts.EpsTarget
        info.converged = true;
        break;
    end

    % Adaptive eps update (Garanzha eq. 6 simplified):
    %   if mesh valid: aggressively reduce eps
    %   if mesh invalid: reduce more cautiously
    if minDJ > 0
        eps_k = max(opts.EpsTarget, eps_k * 0.3);
    else
        eps_k = max(opts.EpsTarget, eps_k * 0.7);
    end
end

msh_out = msh_ref;
msh_out.p = p_phys;
end


% =========================================================================
function [F, grad_int] = energy_and_grad(x, p_phys, int_v, t, Z, volref, lambda, eps_)
nint  = numel(int_v);
p     = p_phys;
p(:, int_v) = reshape(x, 2, nint);
nt    = size(t, 2);

F        = 0;
grad_p   = zeros(size(p));

for it = 1:nt
    v   = t(:, it);
    U   = p(:, v);              % 2 x 3
    Z_t = Z(:, :, it);          % 3 x 2
    J   = U * Z_t;              % 2 x 2

    a = J(1,1); b = J(1,2);
    c = J(2,1); d = J(2,2);

    detJ  = a*d - b*c;
    trJTJ = a*a + b*b + c*c + d*d;

    sqrtED  = sqrt(eps_^2 + detJ^2);
    chi     = (detJ + sqrtED) / 2;
    chi_pr  = (1 + detJ / sqrtED) / 2;     % d chi / d D

    % f_eps = trJTJ / chi
    % g_eps = (detJ^2 + 1) / chi
    f_val = trJTJ / chi;
    g_val = (detJ^2 + 1) / chi;
    F     = F + (f_val + lambda * g_val) * volref(it);

    % Gradients w.r.t. J entries (a,b,c,d)
    %   df/dJ = 2 J / chi  -  (trJTJ / chi^2) * chi_pr * d(detJ)/dJ
    %   dg/dJ = 2 detJ / chi * d(detJ)/dJ  -  ((detJ^2+1)/chi^2)*chi_pr*d(detJ)/dJ
    %   d(detJ)/dJ = [ d  -c ; -b  a ]
    chi2     = chi * chi;
    dD_dJ    = [ d, -c; -b, a ];

    df_dJ = (2 * J) / chi - (trJTJ / chi2) * chi_pr * dD_dJ;
    dg_dJ = (2 * detJ / chi) * dD_dJ - ((detJ^2 + 1) / chi2) * chi_pr * dD_dJ;

    dE_dJ = (df_dJ + lambda * dg_dJ) * volref(it);

    % Chain rule J = U * Z_t  -->  dE/dU = dE/dJ * Z_t^T
    dE_dU = dE_dJ * Z_t';       % 2 x 3

    for j = 1:3
        grad_p(:, v(j)) = grad_p(:, v(j)) + dE_dU(:, j);
    end
end

grad_int = reshape(grad_p(:, int_v), [], 1);
end


% =========================================================================
function [minDJ, allDJ] = min_det_J(p, t, Z)
nt = size(t, 2);
allDJ = zeros(nt, 1);
for it = 1:nt
    v = t(:, it);
    U = p(:, v);
    J = U * Z(:, :, it);
    allDJ(it) = J(1,1)*J(2,2) - J(1,2)*J(2,1);
end
minDJ = min(allDJ);
end


% =========================================================================
function [x, F] = bfgs_minimize(fun, x0, max_iter, gtol, verbose)
% Hand-coded BFGS with Armijo backtracking (same template as optimize_shape).
n = numel(x0);
x = x0(:);
H = eye(n);

[F, g] = fun(x);
if verbose
    fprintf('    inner iter %4d  F=%.4e |g|=%.3e\n', 0, F, norm(g));
end

for iter = 1:max_iter
    if norm(g) < gtol, break; end
    d = -H * g;
    dg = d' * g;
    if dg >= 0
        H = eye(n);
        d = -g;
        dg = d' * g;
    end
    alpha = 1;
    success = false;
    Fn = inf; gn = g;
    for bt = 1:60
        xn = x + alpha * d;
        [Fn, gn] = fun(xn);
        if isfinite(Fn) && Fn <= F + 1e-4 * alpha * dg
            success = true;
            break;
        end
        alpha = alpha / 2;
    end
    if ~success, break; end
    s = xn - x;  y = gn - g;
    sy = s' * y;
    if sy > 1e-12
        rho = 1 / sy;
        I = eye(n);
        H = (I - rho*s*y') * H * (I - rho*y*s') + rho * (s*s');
    end
    x  = xn;
    F  = Fn;
    g  = gn;
end
end


% =========================================================================
function is_bnd = boundary_mask(msh)
np = size(msh.p, 2);
t  = double(msh.t) + 1;
t2t = double(msh.t2t);
is_bnd = false(np, 1);
for it = 1:size(t, 2)
    for j = 1:size(t2t, 1)
        if t2t(j, it) < 0
            v = t(:, it);
            v(j) = [];
            is_bnd(v) = true;
        end
    end
end
end
