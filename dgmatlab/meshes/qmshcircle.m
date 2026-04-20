function msh=qmshcircle(n)

if nargin<1, n=1; end

s=(0:n)'/n;
e=ones(size(s));
a=(2+pi/8)/(4+sqrt(2)); % All coarsest quads same average (curved) side length

ix=reshape(1:(n+1)^2,[n+1,n+1]);
q0=[ix(1:n,1:n),ix(2:n+1,1:n),ix(1:n,2:n+1),ix(2:n+1,2:n+1)];
q0=reshape(q0,[],4);

phi=linspace(0,pi/4,n+1)';
x1=a*e; y1=a*s;
x2=cos(phi); y2=sin(phi);
x3=linspace(a,1,n+1)'; y3=0*x3;
x4=linspace(a,1/sqrt(2),n+1)'; y4=x4;

X1=(1-s)*x1'+s*x2';
Y1=(1-s)*y1'+s*y2';
X2=x3*(1-s')+x4*s';
Y2=y3*(1-s')+y4*s';
X12=(1-s)*(1-s')*X1(1,1)+s*(1-s')*X1(end,1)+(1-s)*s'*X1(1,end)+s*s'*X1(end,end);
Y12=(1-s)*(1-s')*Y1(1,1)+s*(1-s')*Y1(end,1)+(1-s)*s'*Y1(1,end)+s*s'*Y1(end,end);
X=X1+X2-X12;
Y=Y1+Y2-Y12;
p=[X(:),Y(:)];

a1=linspace(0,a,n+1)';
[xx,yy]=ndgrid(a1,a1);
p0=[xx(:),yy(:)];

p=[p; p(:,[2,1]); p0];
q=[q0; q0(:,[1,3,2,4])+(n+1)^2; q0+2*(n+1)^2];

p1=p;
p2=[-p(:,2),p(:,1)];
p3=-p;
p4=-p2;

p=[p1;p2;p3;p4];
q0=q;
for i=1:3
    q=[q; q0+3*(n+1)^2*i];
end
[p,q]=fixmesh(p,q,1e-8);

bndexpr={'true'};
fd=@(p) sqrt(sum(p.^2,2))-1;
msh=ml2msh(p,q,bndexpr,fd);
msh=mshcurved(msh,1);
