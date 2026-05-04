function [msh_final, info] = winslow_homotopy(msh_ref, msh_target, n_steps, ws_opts)
% WINSLOW_HOMOTOPY  Beta-homotopy continuation for elliptic Winslow smoothing.
%
% Steps the boundary condition gradually from the reference mesh's boundary
% to the target mesh's boundary in n_steps stages, running damped-Picard
% Winslow at each stage with the previous stage's result as the initial
% guess. This is the technique described in Fortunato's thesis eq. (2.26):
% when standard Picard or Newton fails to converge directly on a large
% boundary deformation, a homotopy in the boundary condition produces a
% sequence of well-behaved sub-problems.
%
% Inputs:
%   msh_ref     reference mesh whose boundary is the homotopy "start"
%               (typically the shape-opt'd Tutte output). Its p1 also
%               provides the initial guess for the first stage.
%   msh_target  mesh whose boundary is the homotopy "end" (typically the
%               physical mesh whose boundary we ultimately want).
%   n_steps     number of homotopy stages. Each stage advances beta by
%               1/n_steps. (10-20 is typical.)
%   ws_opts     options struct passed to elliptic_smoothing_damped at each
%               stage (alpha, maxiter, tol). Reasonable defaults applied.
%
% Output:
%   msh_final   mesh after the final stage (beta = 1).
%   info        struct: .converged_steps, .residuals (per stage),
%                       .min_J_per_step, .completed (bool).

if nargin < 3, n_steps = 10; end
if nargin < 4, ws_opts = struct(); end
if ~isfield(ws_opts,'alpha'),   ws_opts.alpha   = 0.5; end
if ~isfield(ws_opts,'maxiter'), ws_opts.maxiter = 100; end
if ~isfield(ws_opts,'tol'),     ws_opts.tol     = 1e-5; end

is_bnd_dg   = boundary_dg_mask(msh_ref);
ref_bnd_p1  = msh_ref.p1;
phys_bnd_p1 = msh_target.p1;

current_msh = msh_ref;
info.converged_steps = 0;
info.completed       = false;
info.min_J_per_step  = nan(n_steps, 1);

for step = 1:n_steps
    beta = step / n_steps;
    fprintf('\n  === homotopy step %d/%d  beta=%.3f ===\n', step, n_steps, beta);

    target_p1 = current_msh.p1;
    target_p1(is_bnd_dg) = (1 - beta) * ref_bnd_p1(is_bnd_dg) + ...
                                beta  * phys_bnd_p1(is_bnd_dg);

    try
        current_msh = elliptic_smoothing_damped(current_msh, target_p1, ...
            false, [], ws_opts);
        info.converged_steps = step;
        J = geojac(current_msh);
        info.min_J_per_step(step) = min(J(:));
        if min(J(:)) < 0
            fprintf('  WARNING: step %d produced negative Jacobian (%.3e). Continuing.\n', ...
                step, min(J(:)));
        end
    catch ME
        fprintf('  step %d FAILED hard: %s\n', step, ME.message);
        break;
    end
end

if info.converged_steps == n_steps
    info.completed = true;
end
msh_final = current_msh;
end


% =========================================================================
function m = boundary_dg_mask(msh)
data = dginit(msh);
[~, edg] = getbndnodes(msh, data);
[ns, d, nt] = size(msh.p1);
m2 = false(ns, nt);
m2(edg) = true;
m = repmat(reshape(m2, ns, 1, nt), 1, d, 1);
end
