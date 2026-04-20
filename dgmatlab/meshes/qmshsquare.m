function msh=qmshsquare(m,n,periodic)

if nargin<1, m=1; end
if nargin<2, n=m; end
if nargin<3 | isempty(periodic), periodic=false; end
if prod(size(periodic))==1, periodic=periodic([1,1]); end

[x,y]=ndgrid((0:m)/m,(0:n)/n);
p=[x(:),y(:)];
ix=reshape(1:prod(size(x)),size(x));
q=[ix(1:end-1,1:end-1),ix(2:end,1:end-1),ix(1:end-1,2:end),ix(2:end,2:end)];
q=reshape(q,[],4);

[t2t,t2n]=mkt2t(q);

bndexpr={'all(p(:,2)<1e-3)','all(p(:,1)>1-1e-3)', ...
         'all(p(:,2)>1-1e-3)','all(p(:,1)<1e-3)'};
periodic_map={2,'p(:,2)',4,'p(:,2)';
              1,'p(:,1)',3,'p(:,1)'};

msh=ml2msh(p,q,bndexpr,[],[],periodic_map(periodic,:));
