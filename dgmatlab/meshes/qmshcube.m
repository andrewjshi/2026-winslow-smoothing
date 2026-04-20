function msh=qmshcube(m,n,o,periodic)

if nargin<1, m=1; end
if nargin<2, n=m; end
if nargin<3, o=n; end
if nargin<4 | isempty(periodic), periodic=false; end
if prod(size(periodic))==1, periodic=periodic([1,1,1]); end

[x,y,z]=ndgrid((0:m)/m,(0:n)/n,(0:o)/o);
p=[x(:),y(:),z(:)];
ix=reshape(1:prod(size(x)),size(x));
q=cat(3,ix(1:m,1:n,1:o),ix(2:m+1,1:n,1:o),ix(1:m,2:n+1,1:o),ix(2:m+1,2:n+1,1:o), ...
   ix(1:m,1:n,2:o+1),ix(2:m+1,1:n,2:o+1),ix(1:m,2:n+1,2:o+1),ix(2:m+1,2:n+1,2:o+1));
q=reshape(q,[],8);

bndexpr={'p(:,1)<1e-3','p(:,1)>1-1e-3', ...
         'p(:,2)<1e-3','p(:,2)>1-1e-3', ...
         'p(:,3)<1e-3','p(:,3)>1-1e-3'};
periodic_map={1,'p(:,[2,3])',2,'p(:,[2,3])';
              3,'p(:,[1,3])',4,'p(:,[1,3])';
              5,'p(:,[1,2])',6,'p(:,[1,2])'};

msh=ml2msh(p,q,bndexpr,[],[],periodic_map(periodic,:));
