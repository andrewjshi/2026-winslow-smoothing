function msh=mshnaca2dm

% Requires distmesh

hs=[0.01,0.025,0.025,1];

a=.12/.2*[0.2969,-0.1260,-0.3516,0.2843,-0.1015];
theta=0;
fix=[0,0;1.0089304129,0;-4.75,0;5.25,0];

xx=linspace(fix(1,1),fix(2,1),10);
yy=a(1)*sqrt(xx)+polyval([a(5:-1:2),0],xx);
pnaca=[xx(:),yy(:)];

fd=inline('ddiff(ddiff(dcircle(p,0.25,0,5),p(:,2)),dnaca(p,a,theta))','p','a','theta','hs','pnaca');
fh=inline('min(min(min(pntsrc(p,[0,0],hs(1)),pntsrc(p,[1,0],hs(2))),hs(3)+0.3*abs(dpoly(p,pnaca))),hs(4))','p','a','theta','hs','pnaca');
box=[-4.75,0;5.25,5];

clear pars
pars.popctrl=[200,500];
rand('state',123);
[p,t]=dm2d(fd,fh,min(hs(:)),box,fix,pars,a,theta,hs,pnaca);
p0=p; t0=t;

p=[p0; p0(:,1),-p0(:,2)];
t=[t0;t0+size(p0,1)];
[p,t]=fixmesh(p,t,1e-8);

fd=inline('ddiff(dcircle(p,0.25,0,5),dnaca(p,a,theta))','p','a','theta');
fdargs={a,theta};

bndexpr={'sqrt((p(:,1).^2)+p(:,2).^2)<2','true'};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,[1,2]);
msh=mshfix(msh);
