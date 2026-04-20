function msh=mshbump1(mul)
% BUMP1 (for Euler shock)

if nargin<1, mul=1; end

[p1,e1,t1]=squaremesh(1*mul+1,3*mul+1);
[p2,e2,t2]=squaremesh(2*mul+1,3*mul+1);
[p3,e3,t3]=squaremesh(3*mul+1,3*mul+1);

l1=sqrt(3^2-2.75^2);
f1=3-l1;
f2=3+l1;

p1(:,1)=p1(:,1)*f1;
p1(:,2)=p1(:,2)*5;

p2(:,1)=p2(:,1)*(f2-f1)+f1;
y0=sqrt(3^2-(p2(:,1)-3).^2)-2.75;
p2(:,2)=p2(:,2).*(5-y0)+y0;

p3(:,1)=p3(:,1)*(10-f2)+f2;
p3(:,2)=p3(:,2)*5;

p=[p1;p2;p3];
t=[t1;t2+size(p1,1);t3+size(p1,1)+size(p2,1)];
[p,t]=fixmesh(p,t);

fd=inline('ddiff(drectangle(p,0,10,0,5),dcircle(p,3,-2.75,3))','p');
fdargs={};

bndexpr={'all(p(:,2)<1e-6)','all(p(:,1)>10-1e-6)', ...
         'all(p(:,2)>5-1e-6)','all(p(:,1)<1e-6)','true'};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,5);
msh=mshfix(msh);
