function msh=mshsquare(m,n,periodic,parity)

if nargin<1 | isempty(m), m=10; end
if nargin<2 | isempty(n), n=m; end
if nargin<3 | isempty(periodic), periodic=false; end
if nargin<4, parity=0; end
if prod(size(periodic))==1, periodic=periodic([1,1]); end

[p,e,t]=squaremesh(m,n,parity);
[t2t,t2n]=mkt2t(t);

bndexpr={'all(p(:,2)<1e-3)','all(p(:,1)>1-1e-3)', ...
         'all(p(:,2)>1-1e-3)','all(p(:,1)<1e-3)'};
periodic_map={2,'p(:,2)',4,'p(:,2)';
              1,'p(:,1)',3,'p(:,1)'};

msh=ml2msh(p,t,bndexpr,[],[],periodic_map(periodic,:));
msh=mshcurved(msh,[]);
msh=mshfix(msh);
