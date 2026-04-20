function msh1 = qmshrefine(msh, edge_refine)
% Unstructured high-order quad-mesh refinement

if size(msh.p,1) ~= 2, error('Must be 2-D'); end
if msh.eltype ~= t_block, error('Must be quad mesh'); end

p = msh.p';
q = msh.t' + 1;
q2q = msh.t2t' + 1;
p1 = msh.p1;
s0 = msh.s0;

ns = size(p1,1);
np = size(p,1);
nq = size(q,1);
porder = double(msh.porder);

data = dginit(msh);
egix = double(data.egix + 1);
V0 = plegendre(2*s0 - 1, porder);
Vmid = plegendre(0, porder) / V0;

% List of refinement patterns - edges refined, ref p, ref q
refine_patterns = {
    [0,0,1,1], [-1 -1;1 -1;-1 1;1 1;0 -1;0 1], [5 2 6 4;1 5 3 6];
    [1,1,0,0], [-1 -1;1 -1;-1 1;1 1;-1 0;1 0], [5 6 3 4;1 2 5 6];
    [1,0,1,0], [-1 -1;1 -1;-1 1;1 1;0 -1;-1 0;0 0], [6 7 3 4;1 5 6 7;5 2 7 4];
    [1,0,0,1], [-1 -1;1 -1;-1 1;1 1;-1 0;0 1;0 0], [7 2 6 4;1 2 5 7;5 7 3 6];
    [0,1,1,0], [-1 -1;1 -1;-1 1;1 1;0 -1;1 0;0 0], [7 6 3 4;1 5 3 7;5 2 7 6];
    [0,1,0,1], [-1 -1;1 -1;-1 1;1 1;1 0;0 1;0 0], [7 5 6 4;1 2 7 5;1 7 3 6];
    [1,1,1,1], [-1 -1;1 -1;-1 1;1 1;0 -1;-1 0;1 0;0 1;0 0], [1 5 6 9;5 2 9 7;6 9 3 8;9 7 8 4];
                  };

% Build refinement interpolation operators
Vs = {};
for iref = 1:size(refine_patterns,1)
    msh1 = ml2msh(refine_patterns{iref,2}, refine_patterns{iref,3});
    msh1 = nodealloc(msh1, porder);

    p2 = (porder+1)^2;
    Vs{iref} = zeros(p2,p2,size(msh1.t,2));
    for iq = 1:size(msh1.t,2)
        Vx = plegendre(msh1.p1(:,1,iq), porder) / V0;
        Vy = plegendre(msh1.p1(:,2,iq), porder) / V0;
        Vs{iref}(:,:,iq) = reshape(repmat(Vx,[1,1,porder+1]) .* ...
                                   repmat(permute(Vy,[1,3,2]),[1,porder+1,1]), [p2,p2]);
    end
end

facemap = mkfacemap(2,1);

% Find all edges
e0 = [q(:,facemap(:,1)); q(:,facemap(:,2)); q(:,facemap(:,3)); q(:,facemap(:,4))];
[e,eix,ejx] = unique(sort(e0,2), 'rows');
ne = size(e,1);

% Edges to refine (input + a different representation)
elrefind = edge_refine;
refind = logical(accumarray(ejx, elrefind(:), [ne,1]));
elrefind = reshape(refind(ejx), nq, 4);

% Ensure consistent refinements:
%   - If 3 edges marked, also mark the 4th
%   - If 1 edge marked, also mark the first edge clockwise
% Repeat until no 3's and 1's left
while 1
    i3 = find(sum(elrefind,2) == 3);
    elrefind(i3,:) = true;
    if isempty(i3)
        i1 = find(sum(elrefind,2) == 1);
        if isempty(i1)
            break;
        end
        j = find(elrefind(i1(1),:));
        kmap = [3,4,2,1];
        elrefind(i1(1),kmap(j)) = true;
    end
    % Mark all neighbors of extra refined elements
    refind = logical(accumarray(ejx, elrefind(:), [ne,1]));
    elrefind = reshape(refind(ejx), nq, 4);
end

% Determine # of new elements and whether midpoints are needed
needelmid = logical(zeros(nq,1));
nnew = 0;
for iq = 1:nq
    if sum(elrefind(iq,:)) == 2
        if ismember(elrefind(iq,:), [1,1,0,0; 0,0,1,1], 'rows')
            nnew = nnew + 1;
        elseif ismember(elrefind(iq,:), [1,0,1,0; 1,0,0,1;
                                0,1,1,0; 0,1,0,1], 'rows')
            needelmid(iq) = true;
            nnew = nnew + 2;
        end
    elseif sum(elrefind(iq,:)) == 4
        needelmid(iq) = true;
        nnew = nnew + 3;
    end
end

% Compute new nodes (split edges and compute element midpoints)
pnew = zeros(sum(refind), 2);
[qix,jix] = find(elrefind);
refindix = find(refind);
for ii = 1:numel(refindix)
    ie = refindix(ii);
    iq = mod(eix(ie)-1, nq) + 1;
    j = floor(((eix(ie)-1) / nq)) + 1;
    pnew(ii,:) = Vmid * p1(egix(:,j),:,iq);
end

newmid = zeros(sum(needelmid), 2);
qix = find(needelmid);
for ii = 1:numel(qix)
    iq = qix(ii);
    for dim = 1:2
        newmid(ii,dim) = Vmid * reshape(p1(:,dim,iq), [porder+1,porder+1]) * Vmid';
    end
end

pp = [p; pnew; newmid];

% ixedge = indices of new edge nodes
% ixmid  = indices of new mid nodes
ixedge = zeros(size(refind));
ixedge(refind) = 1:sum(refind);
ixedge = np + reshape(ixedge(ejx), nq, 4);
ixmid(needelmid) = np + sum(refind) + (1:sum(needelmid));

% Finally do all of the refinements
% Replace directly in q + add new elements in qnew
q2qnew = zeros(nnew,4);
qnew = zeros(nnew,4);
qnewix = 0;
p1new = zeros(ns, 2, nnew);
for iq = 1:nq
    cp1 = p1(:,:,iq);
    if sum(elrefind(iq,:)) == 2
        if isequal(elrefind(iq,:), [0,0,1,1]) % split across, 2 quads
            qnewix = qnewix + 1;
            qnew(qnewix,:) = [q(iq,1),ixedge(iq,3),q(iq,3),ixedge(iq,4)];
            q(iq,:) = [ixedge(iq,3),q(iq,2),ixedge(iq,4),q(iq,4)];
            
            p1(:,:,iq) = Vs{1}(:,:,1) * cp1;
            p1new(:,:,qnewix) = Vs{1}(:,:,2) * cp1;

            q2qnew(qnewix,:) = q2q(iq,:);
            q2qnew(qnewix,2) = 1;
            q2q(iq,1) = 1;
        elseif isequal(elrefind(iq,:), [1,1,0,0]) % split across, 2 quads
            qnewix = qnewix + 1;
            qnew(qnewix,:) = [q(iq,1),q(iq,2),ixedge(iq,1),ixedge(iq,2)];
            q(iq,:) = [ixedge(iq,1),ixedge(iq,2),q(iq,3),q(iq,4)];

            p1(:,:,iq) = Vs{2}(:,:,1) * cp1;
            p1new(:,:,qnewix) = Vs{2}(:,:,2) * cp1;
            
            q2qnew(qnewix,:) = q2q(iq,:);
            q2qnew(qnewix,4) = 1;
            q2q(iq,3) = 1;
        elseif isequal(elrefind(iq,:), [1,0,1,0]) % split corner, 3 quads
            qnew(qnewix+(1:2),:) = ...
                [q(iq,1),ixedge(iq,3),ixedge(iq,1),ixmid(iq);
                 ixedge(iq,3),q(iq,2),ixmid(iq),q(iq,4)];
            q(iq,:) = [ixedge(iq,1),ixmid(iq),q(iq,3),q(iq,4)];
            
            p1(:,:,iq) = Vs{3}(:,:,1) * cp1;
            p1new(:,:,qnewix+1) = Vs{3}(:,:,2) * cp1;
            p1new(:,:,qnewix+2) = Vs{3}(:,:,3) * cp1;

            q2qnew(qnewix+1,:) = q2q(iq,:);
            q2qnew(qnewix+2,:) = q2q(iq,:);
            q2qnew(qnewix+1,[2,4]) = 1;
            q2qnew(qnewix+2,[1,4]) = 1;
            q2q(iq,[2,3]) = 1;
            
            qnewix = qnewix + 2;
        elseif isequal(elrefind(iq,:), [1,0,0,1]) % split corner, 3 quads
            qnew(qnewix+(1:2),:) = ...
                [q(iq,1),q(iq,2),ixedge(iq,1),ixmid(iq);
                 ixedge(iq,1),ixmid(iq),q(iq,3),ixedge(iq,4)];
            q(iq,:) = [ixmid(iq),q(iq,2),ixedge(iq,4),q(iq,4)];
            
            p1(:,:,iq) = Vs{4}(:,:,1) * cp1;
            p1new(:,:,qnewix+1) = Vs{4}(:,:,2) * cp1;
            p1new(:,:,qnewix+2) = Vs{4}(:,:,3) * cp1;
            
            q2qnew(qnewix+1,:) = q2q(iq,:);
            q2qnew(qnewix+2,:) = q2q(iq,:);
            q2qnew(qnewix+1,[2,4]) = 1;
            q2qnew(qnewix+2,[2,3]) = 1;
            q2q(iq,[1,3]) = 1;
            
            qnewix = qnewix + 2;
        elseif isequal(elrefind(iq,:), [0,1,1,0]) % split corner, 3 quads
            qnew(qnewix+(1:2),:) = ...
                [q(iq,1),ixedge(iq,3),q(iq,3),ixmid(iq);
                 ixedge(iq,3),q(iq,2),ixmid(iq),ixedge(iq,2)];
            q(iq,:) = [ixmid(iq),ixedge(iq,2),q(iq,3),q(iq,4)];

            p1(:,:,iq) = Vs{5}(:,:,1) * cp1;
            p1new(:,:,qnewix+1) = Vs{5}(:,:,2) * cp1;
            p1new(:,:,qnewix+2) = Vs{5}(:,:,3) * cp1;
            
            q2qnew(qnewix+1,:) = q2q(iq,:);
            q2qnew(qnewix+2,:) = q2q(iq,:);
            q2qnew(qnewix+1,[2,4]) = 1;
            q2qnew(qnewix+2,[1,4]) = 1;
            q2q(iq,[1,3]) = 1;

            qnewix = qnewix + 2;
        elseif isequal(elrefind(iq,:), [0,1,0,1]) % split corner, 3 quads
            qnew(qnewix+(1:2),:) = ...
                [q(iq,1),q(iq,2),ixmid(iq),ixedge(iq,2);
                 q(iq,1),ixmid(iq),q(iq,3),ixedge(iq,4)];
            q(iq,:) = [ixmid(iq),ixedge(iq,2),ixedge(iq,4),q(iq,4)];

            p1(:,:,iq) = Vs{6}(:,:,1) * cp1;
            p1new(:,:,qnewix+1) = Vs{6}(:,:,2) * cp1;
            p1new(:,:,qnewix+2) = Vs{6}(:,:,3) * cp1;
            
            q2qnew(qnewix+1,:) = q2q(iq,:);
            q2qnew(qnewix+2,:) = q2q(iq,:);
            q2qnew(qnewix+1,[1,4]) = 1;
            q2qnew(qnewix+2,[2,3]) = 1;
            q2q(iq,[1,3]) = 1;

            qnewix = qnewix + 2;
        end
    elseif sum(elrefind(iq,:)) == 4 % split uniform, 4 quads
        qnew(qnewix+(1:3),:) = ...
            [ixedge(iq,3),q(iq,2),ixmid(iq),ixedge(iq,2);
             ixedge(iq,1),ixmid(iq),q(iq,3),ixedge(iq,4);
             ixmid(iq),ixedge(iq,2),ixedge(iq,4),q(iq,4)];
        q(iq,:) = [q(iq,1), ixedge(iq,3), ixedge(iq,1), ixmid(iq)];

        p1(:,:,iq) = Vs{7}(:,:,1) * cp1;
        p1new(:,:,qnewix+1) = Vs{7}(:,:,2) * cp1;
        p1new(:,:,qnewix+2) = Vs{7}(:,:,3) * cp1;
        p1new(:,:,qnewix+3) = Vs{7}(:,:,4) * cp1;

        q2qnew(qnewix+1,:) = q2q(iq,:);
        q2qnew(qnewix+2,:) = q2q(iq,:);
        q2qnew(qnewix+3,:) = q2q(iq,:);
        q2qnew(qnewix+1,[1,4]) = 1;
        q2qnew(qnewix+2,[2,3]) = 1;
        q2qnew(qnewix+3,[1,3]) = 1;
        q2q(iq,[2,4]) = 1;
        
        qnewix = qnewix + 3;
    end
end
qq = [q; qnew];
qq2qq = [q2q; q2qnew];
pp1 = cat(3, p1, p1new);

msh1 = ml2msh(pp, qq);
ix = find(qq2qq <= 0);
q2q = msh1.t2t' + 1;
q2q(ix) = qq2qq(ix);
msh1.t2t = int32(q2q' - 1);

msh1 = nodealloc(msh1, porder);
msh1.p1 = pp1;
msh1 = mshcurved(msh1, 'all');

%dgmeshplot_curved(msh1, 2, 1)
