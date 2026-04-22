function WK = compute_WK_simplex(msh, shape)
% Per-element Jacobian W_K of the affine map from an ideal simplex
% to the reference (straight-sided) element, computed from corner
% positions in msh.p only. Returns a (d,d,nt) array.
%
%   shape = 'equilateral' (default)  -- unit equilateral simplex
%   shape = 'right'                  -- unit right simplex (gives W_K
%                                       identical to the dgmatlab
%                                       reference element, so W_K = R
%                                       times constant; useful for
%                                       regression).

if nargin < 2, shape = 'equilateral'; end

d  = size(msh.p, 1);
nt = size(msh.t, 2);

switch shape
case 'equilateral'
    if d == 2
        V = [0, 1, 0.5;
             0, 0, sqrt(3)/2];
    elseif d == 3
        V = [0, 1, 0.5, 0.5;
             0, 0, sqrt(3)/2, sqrt(3)/6;
             0, 0, 0,         sqrt(2/3)];
    else
        error('compute_WK_simplex: dim %d unsupported', d);
    end
case 'right'
    V = [zeros(d,1), eye(d)];
otherwise
    error('compute_WK_simplex: unknown shape ''%s''', shape);
end

M_inv = inv(V(:, 2:end) - V(:, 1));          % (d x d)
t = double(msh.t) + 1;                        % 1-based

WK = zeros(d, d, nt);
for it = 1:nt
    corners = msh.p(:, t(:, it));              % (d x (d+1))
    E = corners(:, 2:end) - corners(:, 1);     % (d x d)
    WK(:, :, it) = E * M_inv;
end
end
