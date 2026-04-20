function sz=mshsize(msh,type)
%MSHSIZE  Compute element sizes in simplex mesh.
%   sz=mshsize(msh,type)
%
%   type = 1 : Average edge lengths at nodes, interpolate linearly (default)
%   type = 2 : Max edge length of element, constant within each element
%   type = 3 : Min edge length of element, constant within each element
%
%   sz is interpolated to full discontinuous basis of degree msh.porder.

if ~isfield(msh,'porder')
  error('Need msh.porder and msh.s for interpolation, run nodealloc first.');
end

if nargin<2, type=1; end

p=msh.p';
t=msh.t'+1;
dim=size(t,2)-1;
np=size(p,1);
nt=size(t,1);

edge=zeros(0,dim);
facemap=mkfacemap(dim);
for i=1:size(facemap,2)
  edge=[edge; t(:,facemap(:,i))];
end

if type==1
  edge=unique(sort(edge,2),'rows');
  L=sqrt(sum((p(edge(:,1),:)-p(edge(:,2),:)).^2,2));
  
  edge=double(edge);
  sz0=full(sparse(edge,1,L(:,[1,1]),np,1))./ ...
      full(sparse(edge,1,1,np,1));
  
  sz1=reshape(sz0(t(:,[end,1:end-1]))',dim+1,1,nt);
  sz=dginterp(sz1,1,msh.porder,dim);
elseif type==2
  L=sqrt(sum((p(edge(:,1),:)-p(edge(:,2),:)).^2,2));
  sz0=max(reshape(L,nt,[]),[],2);
  sz=repmat(permute(sz0',[1,3,2]),[size(msh.s,1),1,1]);
elseif type==3
  L=sqrt(sum((p(edge(:,1),:)-p(edge(:,2),:)).^2,2));
  sz0=min(reshape(L,nt,[]),[],2);
  sz=repmat(permute(sz0',[1,3,2]),[size(msh.s,1),1,1]);  
else
  error('Unknown element size type.');
end
