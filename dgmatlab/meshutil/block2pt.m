function [p,t]=block2pt(X,ind)

if nargin>=2
  XT=block2tets(X,ind);
else
  XT=block2tets(X);
end
dim = size(XT,1);
p0=reshape(XT,dim,[])';
t0=reshape(1:prod(size(XT))/dim,dim+1,size(XT,3))';
[p,t]=fixmesh(p0,t0);
