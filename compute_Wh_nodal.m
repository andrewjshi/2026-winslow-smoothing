function W_h = compute_Wh_nodal(msh, shape, projection, weighting)
% Lumped L2 projection of the per-element constant W_K to a nodal
% field W_h, returned as (ns, d, d, nt) replicated on the DG nodes.
%
%   shape = 'equilateral' (default) | 'right'         -- ideal element
%   projection = 'vhp' (default) | 'p1'
%   weighting per element w(K):
%       'area'     (default) w(K) = |K|, classical lumped L2 (notes Eq. 9)
%       'eta'                w(K) = eta_K  (boosts distorted elements)
%       'etam1'              w(K) = max(eta_K - 1, eps)
%       'eta2'               w(K) = eta_K^2
%
% At each CG node x_n:
%   W_h(x_n) = sum_{K ni x_n} w(K) W_K  /  sum_{K ni x_n} w(K).

if nargin < 2, shape = 'equilateral'; end
if nargin < 3, projection = 'vhp'; end
if nargin < 4, weighting = 'area'; end

d  = size(msh.p, 1);
nt = size(msh.t, 2);
ns = size(msh.s, 1);

W_K = compute_WK_simplex(msh, shape);
t = double(msh.t) + 1;

% Element areas |K| and distortion eta_K (from straight-sided corners)
absK = zeros(1, nt);
etaK = zeros(1, nt);
for it = 1:nt
    corners = msh.p(:, t(:, it));
    E = corners(:, 2:end) - corners(:, 1);
    detE = det(E);
    absK(it) = abs(detE) / factorial(d);
    % eta_K from the 2x2 edge matrix relative to equilateral ideal
    % (same formula as qualmetrics2d but on straight-sided corners).
    Eb = E * [1 -1/sqrt(3); 0 2/sqrt(3)];   % map right-simplex to equilateral-ref
    etaK(it) = sum(Eb(:).^2) / (d * abs(det(Eb))^(2/d));
end

switch weighting
case 'area',  wK = absK;
case 'eta',   wK = etaK;
case 'etam1', wK = max(etaK - 1, eps);
case 'eta2',  wK = etaK.^2;
otherwise,    error('unknown weighting ''%s''', weighting);
end
absK = wK;  % rebind: lumped sums below use absK as the weight

switch projection
case 'vhp'
    % Lumped L2 projection onto the full CG space V_h^p via mapdg2cg.
    q = mapdg2cg(msh.p1);                    % q(i + ns*(it-1)) = CG index
    n_cg = max(q);

    numW = zeros(d, d, n_cg);
    denW = zeros(1, n_cg);
    for it = 1:nt
        w_contrib = absK(it) * W_K(:, :, it);
        for i = 1:ns
            cg = q(i + ns*(it-1));
            numW(:, :, cg) = numW(:, :, cg) + w_contrib;
            denW(cg) = denW(cg) + absK(it);
        end
    end
    W_CG = numW ./ reshape(denW, 1, 1, []);

    W_h = zeros(ns, d, d, nt);
    for it = 1:nt
        for i = 1:ns
            cg = q(i + ns*(it-1));
            W_h(i, :, :, it) = W_CG(:, :, cg);
        end
    end

case 'p1'
    % Lumped L2 projection onto linear vertices only, then P1 interp to DG.
    npv = size(msh.p, 2);
    numW = zeros(d, d, npv);
    denW = zeros(1, npv);
    for it = 1:nt
        w_contrib = absK(it) * W_K(:, :, it);
        for k = 1:size(t, 1)
            vi = t(k, it);
            numW(:, :, vi) = numW(:, :, vi) + w_contrib;
            denW(vi) = denW(vi) + absK(it);
        end
    end
    W_CG = numW ./ reshape(denW, 1, 1, []);

    S = msh.s(:, 1:d);
    phi_lin = [1 - sum(S, 2), S];            % (ns, d+1)

    W_h = zeros(ns, d, d, nt);
    for it = 1:nt
        corners_idx = t(:, it);
        Wc = W_CG(:, :, corners_idx);
        for i = 1:ns
            acc = zeros(d, d);
            for k = 1:(d+1)
                acc = acc + phi_lin(i, k) * Wc(:, :, k);
            end
            W_h(i, :, :, it) = acc;
        end
    end

otherwise
    error('compute_Wh_nodal: unknown projection ''%s''', projection);
end
end
