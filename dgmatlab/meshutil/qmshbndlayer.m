function [msh1,ix_bndel] = qmshbndlayer(msh, bnds, nlayers)
% Unstructured high-order quad-mesh refinement

if nargin >= 3
    msh1 = msh;
    ix_bndel = [];
    for ilayer = 1:nlayers
        [msh1,ix_bndel_new] = qmshbndlayer(msh1, bnds);
        ix_bndel = unique([ix_bndel(:)', ix_bndel_new(:)', size(msh.t,2)+1:size(msh1.t,2)]);
    end
    return;
end

q = msh.t' + 1;
facemap = mkfacemap(2,1);
e0 = [q(:,facemap(:,1)); q(:,facemap(:,2)); q(:,facemap(:,3)); q(:,facemap(:,4))];

% Mark boundary elements for anisotropic boundary layer refinement
[el, j] = find(ismember(-msh.t2t', bnds));
% Mark any edges with one node on the refinement boundary
bndedges = [];
for cj = 1:4
    ix = j == cj;
    bndedges = [bndedges; q(el(ix), facemap(:,cj))];
end
bndnodes = unique(bndedges(:));
elrefind = reshape(sum(ismember(e0, bndnodes), 2) == 1, size(q));

msh1 = qmshrefine(msh, elrefind);
ix_bndel = find(any(elrefind,2));
