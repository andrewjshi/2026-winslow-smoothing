function msh=mshterahedron(n)

if nargin<1, n=2; end

[p,e,t]=tetrahedronmesh(n);
bndexpr={'all(p(:,1)<1e-3)','all(p(:,2)<1e-3)', ...
         'all(p(:,3)<1e-3)','all(sum(p,2)>1-1e-3)'};
msh=ml2msh(p,t,bndexpr);
msh=mshcurved(msh,[]);
