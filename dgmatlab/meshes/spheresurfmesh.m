function [p,t]=spheresurfmesh(n)

if nargin<1, n=1; end

[s,t]=lagrangepnts(n);

theta=1-s(:,2);
phi=s(:,1)./max(theta,eps);
phi(end)=0;

x=cos(pi/2*phi).*sin(pi/2*theta);
y=sin(pi/2*phi).*sin(pi/2*theta);
z=cos(pi/2*theta);

p=[ x, y, z;
   -x, y, z;
    x,-y, z;
   -x,-y, z;
    x, y,-z;
   -x, y,-z;
    x,-y,-z;
   -x,-y,-z];

np=size(x,1);

t1=t;
t2=t; t2(:,[1,2])=t2(:,[2,1]);

t=[t1+0*np; t2+1*np; t2+2*np; t1+3*np;
   t2+4*np; t1+5*np; t1+6*np; t2+7*np];

[p,t]=fixmesh(p,t);
