function msh=mshcircleinrect(mul1,mul2,mul3)

if nargin<1, mul1=1; end
if nargin<2, mul2=1; end
if nargin<3, mul3=1; end

m=4*mul1;
n=12*mul1;
o=4*mul2;

s1=(1:m)/m; z1=0*s1;
s2=(1:n)/n; z2=0*s2;
phi=2*pi*(1:4*o)/(4*o);

pv1=[-4+24*s2,20+z1,20-24*s2,-4+z1;
     -4+z2,-4+8*s1,4+z2,4-8*s1]';
pv2=[cos(phi); sin(phi)]';

% Using polymesh
%[p,t]=polymesh({pv1,pv2},[1,1],[1,0;0,1],[8/m,1.3]);

% Using polytrimesh
hfunc=@(p) min(2*2*pi/(4*o)+0.5*dcircle(p,0,0,1),2/mul3+0.5*ramp(drectangle(p,0,10,-1.5,1.5)));
pv1=pv1([1:end,1],:);
pv2=pv2([1:end,1],:);
[p,t]=polytrimesh({pv1,pv2},[0,0],'pq35YQ',hfunc);

fd=@dcircle;
fdargs={0,0,1};

bndexpr={'all(p(:,2)<-4+1e-6)','all(p(:,1)>20-1e-6)', ...
         'all(p(:,2)>4-1e-6)','all(p(:,1)<-4+1e-6)', ...
         'all(sum(p.^2,2)<2^2)'};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,5);
msh=mshfix(msh);
