function msh_out = tutte_embedding(msh, varargin)
% TUTTE_EMBEDDING  Straight-line Tutte barycentric embedding.
%
%   msh_out = tutte_embedding(msh)
%       Uses the unit square as the target convex polygon; boundary
%       vertices are distributed uniformly by arc length around it.
%
%   msh_out = tutte_embedding(msh, 'TargetCorners', corners)
%       corners is (d x nc), vertices of a convex polygon in CCW
%       order; boundary vertices are distributed uniformly along it.
%
%   msh_out = tutte_embedding(msh, 'PreserveBoundary', true)
%       Uses the input boundary vertex positions as the target
%       boundary positions. Useful when the input already has its
%       boundary on the polygon we want.
%
% Acts on msh.p (linear vertex positions) and msh.t (triangle vertex
% indices, 0-based). High-order nodes (msh.p1) are not touched;
% caller should re-run nodealloc.
%
% Tutte (1963): for a 3-connected planar graph whose outer boundary is
% fixed to a convex polygon and whose interior vertices are placed at
% convex combinations of their neighbours, the straight-line drawing
% is a valid planar embedding --- no inverted or degenerate faces.

opts.TargetCorners    = [0 1 1 0; 0 0 1 1];   % unit square, CCW
opts.TargetShape      = '';                    % '' | 'square' | 'circle'
opts.PreserveBoundary = false;
for k = 1:2:numel(varargin)
    opts.(varargin{k}) = varargin{k+1};
end

d  = size(msh.p, 1);
np = size(msh.p, 2);

% 1. Extract the boundary as a cyclic list of linear-vertex indices,
%    oriented CCW.
bnd_v = extract_boundary_loop(msh);
if signed_area(msh.p(:, bnd_v)) < 0
    bnd_v = flip(bnd_v);
end
nbnd  = numel(bnd_v);

% 2. Determine boundary target positions.
if opts.PreserveBoundary
    target_bnd = msh.p(:, bnd_v);
elseif strcmpi(opts.TargetShape, 'circle')
    theta = 2*pi * (0:nbnd-1)' / nbnd;
    target_bnd = [cos(theta)'; sin(theta)'];
elseif strcmpi(opts.TargetShape, 'square') || isempty(opts.TargetShape)
    corners = [0 1 1 0; 0 0 1 1];
    if ~isempty(opts.TargetCorners) && ~isequal(opts.TargetCorners, corners)
        corners = opts.TargetCorners;
    end
    target_bnd = distribute_around_polygon(corners, nbnd);
else
    error('tutte_embedding: unknown TargetShape ''%s''', opts.TargetShape);
end

% 3. Build the graph Laplacian L = D - A of the linear vertex connectivity.
t = double(msh.t) + 1;
I = reshape(t([1 2 3 1 2 3], :), [], 1);
J = reshape(t([2 3 1 3 1 2], :), [], 1);
A = sparse(I, J, 1, np, np);
A = (A + A') > 0;                      % symmetric, binary
deg = full(sum(A, 2));
L = spdiags(deg, 0, np, np) - double(A);

% 4. Partition into boundary / interior and solve L_ii x_int = -L_ib x_bnd.
is_bnd       = false(np, 1);
is_bnd(bnd_v)= true;
int_v        = find(~is_bnd);

L_ii = L(int_v, int_v);
L_ib = L(int_v, bnd_v);

x_bnd = target_bnd';                   % nbnd x d
x_int = L_ii \ (-L_ib * x_bnd);        % nint x d

% 5. Assemble and return.
p_new              = zeros(d, np);
p_new(:, bnd_v)    = target_bnd;
p_new(:, int_v)    = x_int';

msh_out      = msh;
msh_out.p    = p_new;
end


% -------------------------------------------------------------------------

function bnd_v = extract_boundary_loop(msh)
% Boundary linear-vertex indices in cyclic order.
t   = double(msh.t) + 1;
t2t = double(msh.t2t);
nt  = size(t, 2);

% Collect unordered boundary edges (face j of triangle it with t2t < 0).
edges = [];
for it = 1:nt
    for j = 1:size(t2t, 1)
        if t2t(j, it) < 0
            v = t(:, it);
            v(j) = [];
            edges = [edges; v(:)']; %#ok<AGROW>
        end
    end
end
nedges = size(edges, 1);
if nedges == 0
    error('tutte_embedding: no boundary edges found.');
end

% Walk the edges to form a cyclic loop.
bnd_v      = zeros(nedges, 1);
bnd_v(1)   = edges(1, 1);
curr       = edges(1, 2);
used       = false(nedges, 1);
used(1)    = true;

for k = 2:nedges
    bnd_v(k) = curr;
    found = false;
    for ie = 1:nedges
        if used(ie), continue; end
        if edges(ie, 1) == curr
            curr = edges(ie, 2); used(ie) = true; found = true; break;
        elseif edges(ie, 2) == curr
            curr = edges(ie, 1); used(ie) = true; found = true; break;
        end
    end
    if ~found
        error('tutte_embedding: failed to close the boundary loop.');
    end
end
end


function s = signed_area(p)
% Shoelace signed area of a 2D polygon given as (2 x n) vertex list.
n = size(p, 2);
s = 0;
for k = 1:n
    kn = mod(k, n) + 1;
    s = s + p(1, k) * p(2, kn) - p(1, kn) * p(2, k);
end
s = s / 2;
end


function target = distribute_around_polygon(corners, nbnd)
% Place nbnd points on a convex polygon such that nc of them land
% exactly on the polygon's corners, and the rest are distributed
% uniformly along each edge. Result is a genuine nc-cornered polygon
% (not an nbnd-gon inscribed in it).
d  = size(corners, 1);
nc = size(corners, 2);
assert(nbnd >= nc, 'tutte_embedding: need nbnd >= nc corners.');

% Boundary-vertex index that snaps to corner c (1-indexed).
corner_k = round((0:nc-1) * nbnd / nc) + 1;

target = zeros(d, nbnd);
for c = 1:nc
    k0 = corner_k(c);
    k1 = corner_k(mod(c, nc) + 1);
    if k1 <= k0, k1 = k1 + nbnd; end
    n = k1 - k0;                       % boundary vertices on this edge, incl. start
    p0 = corners(:, c);
    p1 = corners(:, mod(c, nc) + 1);
    for i = 0:n-1
        t = i / n;
        k = mod(k0 - 1 + i, nbnd) + 1;
        target(:, k) = p0 + t * (p1 - p0);
    end
end
end
