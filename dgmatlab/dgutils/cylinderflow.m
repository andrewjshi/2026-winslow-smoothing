function u=cylflow(p1,R,U,M)

pinf=1/1.4/M^2;

x=p1(:,1,:);
y=p1(:,2,:);

r2=x.^2+y.^2;

theta=atan2(y,x);

u_r=U*cos(theta).*(1-R^2./r2);
u_theta=-U*sin(theta).*(1+R^2./r2);

r=0*x+1;
u=u_r.*cos(theta)-u_theta.*sin(theta);
v=u_r.*sin(theta)+u_theta.*cos(theta);

p=pinf+.5-.5*(u.^2+v.^2);
r=(p./pinf).^(1/1.4);
E=p./r/.4+.5*(u.^2+v.^2);

u=cat(2,r,r.*u,r.*v,r.*E);
