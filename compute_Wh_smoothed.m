function W_h = compute_Wh_smoothed(msh, shape, n_iter)
% Heavy-smoothed nodal W_h field. Starts from the lumped P1 projection
% used in compute_Wh_nodal(..., 'p1'), then applies n_iter rounds of
% Jacobi-style Laplacian smoothing on the linear-vertex graph, then
% P1-interpolates to DG nodes.
%
%   shape  = 'equilateral' (default) | 'right'
%   n_iter = number of smoothing iterations (n_iter = 0 reduces to the
%            unsmoothed P1 projection).
%
% Each iteration replaces v[i] with (v[i] + sum_{j in N(i)} v[j]) /
% (1 + |N(i)|), i.e. the mean over {i} and its linear-vertex neighbors.
% The fixed point of this operator is a constant field, so increasing
% n_iter monotonically reduces grad(W_h) toward zero, at the cost of
% per-element targeting precision.

if nargin < 2, shape = 'equilateral'; end
if nargin < 3, n_iter = 0; end

d   = size(msh.p, 1);
nt  = size(msh.t, 2);
ns  = size(msh.s, 1);
npv = size(msh.p, 2);

W_K = compute_WK_simplex(msh, shape);
t   = double(msh.t) + 1;

% Element areas (lumped P1 weight)
absK = zeros(1, nt);
for it = 1:nt
    corners = msh.p(:, t(:, it));
    E = corners(:, 2:end) - corners(:, 1);
    absK(it) = abs(det(E)) / factorial(d);
end

% Lumped L^2 projection onto linear vertices: area-weighted mean of W_K
W_v = zeros(d, d, npv);
denW = zeros(1, npv);
for it = 1:nt
    contrib = absK(it) * W_K(:, :, it);
    for k = 1:size(t, 1)
        vi = t(k, it);
        W_v(:, :, vi) = W_v(:, :, vi) + contrib;
        denW(vi) = denW(vi) + absK(it);
    end
end
W_v = W_v ./ reshape(denW, 1, 1, []);

% Linear-vertex graph adjacency (vertices i, j adjacent iff some
% triangle contains both)
I = reshape(t([1 2 3 1 2 3], :), [], 1);
J = reshape(t([2 3 1 3 1 2], :), [], 1);
adj = sparse(I, J, 1, npv, npv);
adj = (adj + adj') > 0;
deg = max(full(sum(adj, 2)), 1);

% Iterative Laplacian smoothing (mean over self + neighbors)
for k = 1:n_iter
    W_new = zeros(d, d, npv);
    for a = 1:d
        for b = 1:d
            v = squeeze(W_v(a, b, :));
            v_neigh = adj * v;
            W_new(a, b, :) = (v + v_neigh) ./ (1 + deg);
        end
    end
    W_v = W_new;
end

% P1 interpolation back to DG nodes
S = msh.s(:, 1:d);
phi_lin = [1 - sum(S, 2), S];

W_h = zeros(ns, d, d, nt);
for it = 1:nt
    corners_idx = t(:, it);
    Wc = W_v(:, :, corners_idx);
    for i = 1:ns
        acc = zeros(d, d);
        for k = 1:(d+1)
            acc = acc + phi_lin(i, k) * Wc(:, :, k);
        end
        W_h(i, :, :, it) = acc;
    end
end
end
