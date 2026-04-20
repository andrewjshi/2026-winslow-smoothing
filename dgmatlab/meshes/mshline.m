function msh=mshline(n,periodic)

if nargin<1 | isempty(n), n=11; end
if nargin<2 | isempty(periodic), periodic=false; end

p=(0:n-1)'/(n-1);
t=[1:n-1; 2:n]';
if periodic
  periodic_map={1,'true',2,'true'};
else
  periodic_map={};
end
msh=ml2msh(p,t,{'p(:,1)<1e-8','p(:,1)>1-1e-8'},[],[],periodic_map);
msh=mshcurved(msh,[]);
