function msh_phys_new = derive_physical_from_slid_circle(msh_r, msh_phys, porder, ref_v)
% DERIVE_PHYSICAL_FROM_SLID_CIRCLE  Given a shape-opt output `msh_r` whose
% boundary linear vertices sit on the unit circle (some of them slid by
% optimize_shape with SlideBoundary='circle'), derive new physical
% boundary positions on `msh_phys`'s domain by proportional-arc-length
% pullback. Returns a copy of `msh_phys` with its boundary linear vertices
% moved to the pulled-back positions and high-order DG nodes regenerated.
%
% `ref_v` is a linear-vertex index of a PINNED boundary vertex (e.g. an
% outer corner) to use as the arc-length reference. Using a pinned vertex
% is essential: its angle on the circle is unchanged by sliding, so the
% cum_arc reference frame doesn't drift.

bnd_v = extract_bnd_loop(msh_r);
if signed_area_local(msh_r.p(:, bnd_v)) < 0
    bnd_v = flip(bnd_v);
end

i_ref = find(bnd_v == ref_v, 1);
if isempty(i_ref)
    error('derive_physical_from_slid_circle: ref_v not found on boundary.');
end
bnd_v = bnd_v([i_ref:end, 1:i_ref-1]);
nbnd = numel(bnd_v);

% Cumulative arc length on the ORIGINAL physical perimeter, walked along
% bnd_v in cyclic order. (Polyline approximation of the smooth boundary
% if applicable; for a polygonal boundary this is exact.)
edge_lens = vecnorm( ...
    msh_phys.p(:, bnd_v([2:end, 1])) - msh_phys.p(:, bnd_v), 2, 1)';
cum_arc = [0; cumsum(edge_lens)];
P = cum_arc(end);

slid_pos = msh_r.p(:, bnd_v);
th = atan2(slid_pos(2,:), slid_pos(1,:))';
th_rel = mod(th - th(1), 2*pi);
th_rel(1) = 0;

s_new = th_rel * P / (2*pi);

new_bnd_pos = zeros(2, nbnd);
for i = 1:nbnd
    s = s_new(i);
    k = find(cum_arc(1:end-1) <= s & s < cum_arc(2:end), 1);
    if isempty(k), k = nbnd; end
    denom = cum_arc(k+1) - cum_arc(k);
    if denom < 1e-12
        frac = 0;
    else
        frac = (s - cum_arc(k)) / denom;
    end
    p1 = msh_phys.p(:, bnd_v(k));
    p2 = msh_phys.p(:, bnd_v(mod(k, nbnd) + 1));
    new_bnd_pos(:, i) = p1 + frac * (p2 - p1);
end

msh_phys_new = msh_phys;
msh_phys_new.p(:, bnd_v) = new_bnd_pos;
msh_phys_new = nodealloc(msh_phys_new, porder);
end


% =========================================================================
function bnd_v = extract_bnd_loop(msh)
t = double(msh.t) + 1;
t2t = double(msh.t2t);
edges = [];
for it = 1:size(t, 2)
    for j = 1:size(t2t, 1)
        if t2t(j, it) < 0
            v = t(:, it); v(j) = [];
            edges = [edges; v(:)']; %#ok<AGROW>
        end
    end
end
nedges = size(edges, 1);
bnd_v = zeros(nedges, 1);
bnd_v(1) = edges(1, 1);
curr = edges(1, 2);
used = false(nedges, 1);
used(1) = true;
for k = 2:nedges
    bnd_v(k) = curr;
    for ie = 1:nedges
        if used(ie), continue; end
        if edges(ie, 1) == curr
            curr = edges(ie, 2); used(ie) = true; break;
        elseif edges(ie, 2) == curr
            curr = edges(ie, 1); used(ie) = true; break;
        end
    end
end
end


function s = signed_area_local(p)
n = size(p, 2);
s = 0;
for k = 1:n
    kn = mod(k, n) + 1;
    s = s + p(1, k) * p(2, kn) - p(1, kn) * p(2, k);
end
s = s / 2;
end
