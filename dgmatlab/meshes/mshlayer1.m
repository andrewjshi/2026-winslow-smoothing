function msh=mshlayer1(m,n)

if nargin<1 | isempty(m), m=10; end
if nargin<2 | isempty(n), n=m; end

[p,e,t]=squaremesh(m,n);
p(:,1)=p(:,1).*p(:,2)+(1-p(:,2)).*(exp(alpha*p(:,1))-1)/(exp(alpha*1)-1);

bndexpr={'all(p(:,2)<1e-3)','all(p(:,1)>1-1e-3)', ...
         'all(p(:,2)>1-1e-3)','all(p(:,1)<1e-3)'};
msh=ml2msh(p,t,bndexpr);
msh=mshcurved(msh,[]);
