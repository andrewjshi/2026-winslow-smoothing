function [CL,CD]=liftdrag(u,msh,data,phys,bndnbr,fcn)

if nargin<6, fcn=@dgnavierstokes; end

cmd=sprintf('postint %d %d',1,bndnbr);
uflux=feval(fcn,u,msh,data,phys,cmd);

dim=size(msh.p,1);
uinf=phys.uinf(2:dim+1);

uinf_norm=uinf/norm(uinf);
D=dot(uflux,uinf_norm);
L=norm(uflux-D*uinf_norm(:));

CD=2*D/norm(uinf)/phys.uinf(1);
CL=2*L/norm(uinf)/phys.uinf(1);
