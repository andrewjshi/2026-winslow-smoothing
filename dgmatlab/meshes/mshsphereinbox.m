function msh=mshsphereinbox(mul1,mul2)

if nargin<1, mul1=2; end
if nargin<2, mul2=2; end

m=2*mul1;
n=3*mul1;
o=1*mul2;

[p1,foo,t1]=squaremesh(m,m);
[p2,foo,t2]=squaremesh(n,m);
t1flip=t1(:,[2,1,3]);
t2flip=t2(:,[2,1,3]);

x1=p1(:,1); y1=p1(:,2); z1=0*x1;
x2=p2(:,1); y2=p2(:,2); z2=0*x2;

p3=[-4+z1, -4+8*x1, -4+8*y1];
p4=[8+z1, -4+8*x1, -4+8*y1];
p5=[-4+12*x2, -4+8*y2, 4+z2];
p6=[-4+12*x2, -4+8*y2,-4+z2];
p7=[-4+12*x2, 4+z2, -4+8*y2];
p8=[-4+12*x2,-4+z2, -4+8*y2];

np1=size(p1,1); np2=size(p2,1);
pp1=[p3;p4;p5;p6;p7;p8];
tt1=[t1flip; t1+np1; t2+2*np1; t2flip+2*np1+np2; t2flip+2*np1+2*np2; t2+2*np1+3*np2];

[p,t]=fixmesh(pp1,tt1);

[p1,t1]=spheresurfmesh(o);
t1=t1(:,[2,1,3]);

np=size(p,1);
p=[p;p1];
t=[t;t1+np];

[p,t]=mesh3D(p,t,[8/4/mul1,1.2]);

fd=@dsphere;
fdargs={0,0,0,1};

bndexpr={'sum(p.^2,2)<2^2','p(:,1)<-4+1e-3','p(:,1)>8-1e-3', ...
         'p(:,2)<-4+1e-3','p(:,2)>4-1e-3','p(:,3)<-4+1e-3','p(:,3)>4-1e-3'};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,1);
msh=mshfix(msh);
