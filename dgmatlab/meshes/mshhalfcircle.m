function msh=mshhalfcircle(mul,hmax)

if nargin<1, mul=2; end
if nargin<2, hmax=5; end

n=4*mul;

phi=pi*(n:-1:0)/n;
pv1=[cos(phi); sin(phi)]';

r=10;
pv2=[r,0; r,r; -r,r; -r,0];

pv=[pv1;pv2];

%[p,t]=polymesh({pv},[1],[1,0],[hmax,1.3]);
[p,t]=polytrimesh({pv([1:end,1],:)},[],'pq25Qa01');

fd=@dcircle;
fdargs={0,0,1};

bndexpr={'all(sum(p.^2,2)<1+1e-6)','all(p(:,1)>10-1e-6)', ...
         'all(p(:,2)>10-1e-6)','all(p(:,1)<-10+1e-6)', ...
         'all(p(:,2)<1e-6)'};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,1);
msh=mshfix(msh);
