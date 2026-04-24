function [msh_out, info] = optimize_shape(msh, varargin)
% OPTIMIZE_SHAPE  Global shape optimisation of the straight-sided
% linear triangulation in msh. Minimises sum_K eta_K^2, where eta_K is
% the Liu-Joe inverse mean-ratio distortion of element K measured
% against an ideal equilateral triangle. Boundary linear vertices are
% pinned; interior linear vertices are the only degrees of freedom.
%
%   msh_out = optimize_shape(msh)
%   msh_out = optimize_shape(msh, 'Display', 'iter')
%   msh_out = optimize_shape(msh, 'MaxIter', 500)
%   [msh_out, info] = optimize_shape(...)
%
% Solver: hand-coded BFGS with Armijo backtracking. Analytic gradient.

display_mode = 'none';
max_iter     = 500;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'display',       display_mode = varargin{k+1};
        case 'maxiterations', max_iter     = varargin{k+1};
        otherwise
            error('optimize_shape: unknown option ''%s''', varargin{k});
    end
end

p = double(msh.p);
t = double(msh.t) + 1;

is_bnd = boundary_mask(msh);
int_v  = find(~is_bnd);
nint   = numel(int_v);

% Equilateral-ideal constants, M_eq = [[1, 1/2], [0, sqrt(3)/2]]
% Q = (a^2 + c^2)*A1 + 2*(a*b + c*d)*A2 + (b^2 + d^2)*A3
A1 = 4/3;
A2 = -2/3;
A3 = 4/3;
det_Minv = 2/sqrt(3);

fun = @(x) eta2_sum_and_grad(x, p, t, int_v, A1, A2, A3, det_Minv);

x0 = reshape(p(:, int_v), [], 1);
[F0, ~] = fun(x0);

[x_opt, F_final, iters] = bfgs_minimize(fun, x0, ...
    'MaxIter', max_iter, 'Display', display_mode);

p(:, int_v) = reshape(x_opt, 2, nint);
msh_out = msh;
msh_out.p = p;

info.F_initial  = F0;
info.F_final    = F_final;
info.iterations = iters;
info.nint       = nint;
end


% ========================================================================
function [x, F, iter] = bfgs_minimize(fun, x0, varargin)
% Simple BFGS with Armijo backtracking line search. fun(x) returns
% [F, gradF]. Inverse-Hessian approximation held dense (OK for small n).

opts.MaxIter = 500;
opts.GradTol = 1e-8;
opts.StepTol = 1e-12;
opts.Display = 'none';
opts.c1      = 1e-4;
for k = 1:2:numel(varargin)
    opts.(varargin{k}) = varargin{k+1};
end

x = x0(:);
n = numel(x);
H = eye(n);        % inverse-Hessian approximation

[F, g] = fun(x);

verbose = strcmpi(opts.Display, 'iter');
if verbose
    fprintf('  iter %4d  F = %.6e  |g| = %.3e\n', 0, F, norm(g));
end

iter = 0;
for iter = 1:opts.MaxIter
    if norm(g) < opts.GradTol
        break;
    end

    d = -H * g;
    dg = d' * g;
    if dg >= 0
        % Not a descent direction; reset to steepest descent.
        H  = eye(n);
        d  = -g;
        dg = d' * g;
    end

    % Armijo backtracking
    alpha = 1;
    Fn = Inf; gn = g;
    success = false;
    for bt = 1:60
        xn       = x + alpha * d;
        [Fn, gn] = fun(xn);
        if isfinite(Fn) && Fn <= F + opts.c1 * alpha * dg
            success = true;
            break;
        end
        alpha = alpha / 2;
    end
    if ~success
        if verbose, fprintf('  line search failed at iter %d\n', iter); end
        break;
    end
    if alpha * norm(d) < opts.StepTol
        break;
    end

    % BFGS update
    s  = xn - x;
    y  = gn - g;
    sy = s' * y;
    if sy > 1e-12
        rho = 1 / sy;
        I = eye(n);
        H = (I - rho * s * y') * H * (I - rho * y * s') + rho * (s * s');
    end

    x = xn;
    F = Fn;
    g = gn;

    if verbose
        fprintf('  iter %4d  F = %.6e  |g| = %.3e  alpha = %.2e\n', ...
                iter, F, norm(g), alpha);
    end
end
end


% ========================================================================
function [F, gradF] = eta2_sum_and_grad(x, p, t, int_v, A1, A2, A3, det_Minv)
% Sum_K eta_K^2 and its gradient w.r.t. the interior vertex positions.
p(:, int_v) = reshape(x, 2, numel(int_v));

nt     = size(t, 2);
grad_p = zeros(size(p));
F      = 0;

for it = 1:nt
    v  = t(:, it);
    p1 = p(:, v(1));
    p2 = p(:, v(2));
    p3 = p(:, v(3));

    a = p2(1) - p1(1);
    b = p3(1) - p1(1);
    c = p2(2) - p1(2);
    d = p3(2) - p1(2);

    detE = a * d - b * c;
    if detE <= 0
        F = Inf;
        gradF = zeros(numel(x), 1);
        return;
    end

    Q   = (a^2 + c^2) * A1 + 2 * (a*b + c*d) * A2 + (b^2 + d^2) * A3;
    D   = 2 * det_Minv * detE;
    eta = Q / D;
    F   = F + eta^2;

    dQ_da = 2 * (a*A1 + b*A2);
    dQ_db = 2 * (a*A2 + b*A3);
    dQ_dc = 2 * (c*A1 + d*A2);
    dQ_dd = 2 * (c*A2 + d*A3);

    dD_da =  2 * det_Minv * d;
    dD_db = -2 * det_Minv * c;
    dD_dc = -2 * det_Minv * b;
    dD_dd =  2 * det_Minv * a;

    invD     = 1 / D;
    deta_da  = (dQ_da - eta * dD_da) * invD;
    deta_db  = (dQ_db - eta * dD_db) * invD;
    deta_dc  = (dQ_dc - eta * dD_dc) * invD;
    deta_dd  = (dQ_dd - eta * dD_dd) * invD;

    coeff = 2 * eta;    % d(eta^2)/d(eta)
    dF_da = coeff * deta_da;
    dF_db = coeff * deta_db;
    dF_dc = coeff * deta_dc;
    dF_dd = coeff * deta_dd;

    grad_p(1, v(1)) = grad_p(1, v(1)) - dF_da - dF_db;
    grad_p(1, v(2)) = grad_p(1, v(2)) + dF_da;
    grad_p(1, v(3)) = grad_p(1, v(3)) + dF_db;
    grad_p(2, v(1)) = grad_p(2, v(1)) - dF_dc - dF_dd;
    grad_p(2, v(2)) = grad_p(2, v(2)) + dF_dc;
    grad_p(2, v(3)) = grad_p(2, v(3)) + dF_dd;
end

gradF = reshape(grad_p(:, int_v), [], 1);
end


% ========================================================================
function is_bnd = boundary_mask(msh)
np  = size(msh.p, 2);
t   = double(msh.t) + 1;
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
