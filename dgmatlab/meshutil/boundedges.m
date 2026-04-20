function e=boundedges(p,t,eltype)
%BOUNDEDGES Find boundary edges of mesh
%   E=BOUNDEDGES(P,T)

if nargin<3, eltype=t_simplex; end

% Form all edges, non-duplicates are boundary edges
edges=mkfaces(t,eltype);
edges=sort(edges,2);
[foo,ix,jx]=unique(edges,'rows');
vec=histc(jx,1:max(jx));
qx=find(vec==1);
e=edges(ix(qx),:);

% if eltype==t_simplex
%     node3=[t(:,3);t(:,2);t(:,1)];
%     node3=node3(ix(qx));
    
%     % Orientation
%     v1=p(e(:,2),:)-p(e(:,1),:);
%     v2=p(node3,:)-p(e(:,1),:);
%     ix=find(v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)>0);
%     e(ix,[1,2])=e(ix,[2,1]);
% end
