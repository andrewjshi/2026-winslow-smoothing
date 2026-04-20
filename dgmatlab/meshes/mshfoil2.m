function msh=mshfoil2(res,nbndlayer,ratiobndlayer)

% res = [x1,x2,x3,y1,ybnd]

if nargin<1, res=[2,2,2,2,2]; end
if nargin<2, nbndlayer=11; end
if nargin<3, ratiobndlayer=5; end

load foil1bnd
p1=p;
np1=size(p1,1);
pp1=interp1(1:np1,p1,1:0.25:np1,'cubic');

fd=inline('dpoly(p,pp)','p','pp');
fdargs={pp1};

[p,t]=foilmesh2(res,nbndlayer,ratiobndlayer);

bndexpr={'all((p(:,1)-.25).^2+p(:,2).^2<1)','true'};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,[]);
msh.tcurved(:)=true; msh.ecurved(:)=true; msh.ncurved(:)=true;
msh=mshfix(msh);
