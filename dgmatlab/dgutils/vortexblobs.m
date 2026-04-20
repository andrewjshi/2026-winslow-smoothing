function u=vortexblobs(p,U0,r0)

gamma=1.4;
K=exp(1.25)*U0/r0;
a=-1.25/r0^2;

[u1,v1]=vortexvelo(p,[0,1],K,a);
[u2,v2]=vortexvelo(p,[0,-1],K,a);

u=u1+u2;
v=v1+v2;

p0=1/gamma;
% To do...

function [u,v]=vortexvelo(p,p0,K,a)

x=p(:,1,:);
y=p(:,2,:);

r=sqrt((x-p0(1)).^2+(y-p0(2)).^2);
th=atan2(y,x);

ix=find(r==0);
r(ix)=1;
utheta=K/2/a*(1-exp(-a*r.^2))./r;
utheta(ix)=0;

u=sin(th).*utheta;
v=cos(th).*utheta;
