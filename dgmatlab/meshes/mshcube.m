function msh=mshcube(m,n,o,periodic)

if nargin<1 | isempty(m), m=6; end
if nargin<2 | isempty(n), n=m; end
if nargin<3 | isempty(o), o=n; end
if nargin<4 | isempty(periodic), periodic=false; end
if prod(size(periodic))==1, periodic=periodic([1,1,1]); end

[p,e,t]=cubemesh(m,n,o);

bndexpr={'all(p(:,1)<1e-3)','all(p(:,1)>1-1e-3)', ...
         'all(p(:,2)<1e-3)','all(p(:,2)>1-1e-3)', ...
         'all(p(:,3)<1e-3)','all(p(:,3)>1-1e-3)'};
periodic_map={1,'p(:,[2,3])',2,'p(:,[2,3])';
              3,'p(:,[1,3])',4,'p(:,[1,3])';
              5,'p(:,[1,2])',6,'p(:,[1,2])'};

msh=ml2msh(p,t,bndexpr,[],[],periodic_map(periodic,:));
msh=mshcurved(msh,[]);
