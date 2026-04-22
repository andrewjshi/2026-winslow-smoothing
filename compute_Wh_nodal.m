function W_h = compute_Wh_nodal(msh, shape, projection, method, weighting)
% L2 projection of the per-element constant W_K to a nodal field
% W_h, returned as (ns, d, d, nt) replicated on the DG nodes.
%
%   shape = 'equilateral' (default) | 'right'         -- ideal element
%   projection = 'vhp' (default) | 'p1'               -- target space
%   method = 'full' (default) | 'lumped'              -- mass matrix
%       'full'   solves M W_h = b with the consistent CG mass matrix M,
%                matching the alpha-projection of the original FP code.
%       'lumped' replaces M by its row-sum diagonal, reducing the solve
%                to weighted nodal averaging (Gemini's Eq. 9 form).
%   weighting per element w(K) (lumped method only):
%       'area'     (default) w(K) = |K|
%       'eta'                w(K) = eta_K
%       'etam1'              w(K) = max(eta_K - 1, eps)
%       'eta2'               w(K) = eta_K^2

if nargin < 2, shape = 'equilateral'; end
if nargin < 3, projection = 'vhp'; end
if nargin < 4, method = 'full'; end
if nargin < 5, weighting = 'area'; end

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
    q = mapdg2cg(msh.p1);                    % q(i + ns*(it-1)) = CG index
    n_cg = max(q);

    switch method
    case 'full'
        % Full (consistent) L2 projection onto V_h^p:  M * W_CG = b,
        % component-wise in W. Assemble M with dgmass_ml and b by
        % quadrature of (phi_i * W_K) over each element.
        data = dginit(msh);
        M = dgmass_ml(msh, data);
        phi_g = permute(data.gfs(:,1,:), [1,3,2]);   % (ns, ngauss)
        J = geojac(msh, data);                        % (ngauss, nt)

        intphi = zeros(ns, nt);
        for it = 1:nt
            mul = data.gw .* J(:,it) / 2;
            intphi(:, it) = phi_g * mul;              % int_K phi_i dx
        end

        W_CG = zeros(d, d, n_cg);
        for a = 1:d
            for b = 1:d
                rhs = zeros(n_cg, 1);
                for it = 1:nt
                    wab = W_K(a, b, it);
                    for i = 1:ns
                        cg = q(i + ns*(it-1));
                        rhs(cg) = rhs(cg) + intphi(i, it) * wab;
                    end
                end
                W_CG(a, b, :) = M \ rhs;
            end
        end

    case 'lumped'
        % Lumped L2 projection: weighted nodal averaging by absK(it).
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

    otherwise
        error('compute_Wh_nodal: unknown method ''%s''', method);
    end

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
