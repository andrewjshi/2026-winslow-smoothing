function msh=mshtriangle(nref)

if nargin<1, nref=3; end

p=[0,0;1,0;0,1];
t=int32([1,2,3]);
[p,t]=uniref(p,t,nref);

bndexpr={'all(p(:,2)<1e-3)','all(p(:,1)+p(:,2)>1-1e-3)', ...
         'all(p(:,1)<1e-3)'};
msh=ml2msh(p,t,bndexpr);
msh=mshcurved(msh,[]);
