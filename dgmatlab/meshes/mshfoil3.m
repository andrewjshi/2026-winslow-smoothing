function msh=mshfoil3(R,res)

% todo: fix hanging nodes

if nargin<1, R=10; end

% res = [x1,x2,x3,y1,ybnd]

if nargin<1, res=[4,4,4,2,2]; end
if nargin<2, nbndlayer=11; end
if nargin<3, ratiobndlayer=5; end

load foil1bnd
p1=p;
np1=size(p1,1);
pp1=interp1(1:np1,p1,1:0.25:np1,'cubic');

fd=inline('dpoly(p,pp)','p','pp');
fdargs={pp1};

[p,t]=foilmesh2(res,nbndlayer,ratiobndlayer,1);

n=16;
phi=2*pi*(1:n)'/n;
pc=R*[cos(phi),sin(phi)];

e=segcollect(boundedges(p,t));

[pp,tt]=polymesh({p(e{1}(1:end-1),:),pc},[1,1],[1,0;1,0],[10,1.9]);
%for ii=1:5
%  p=bndproj(p,t,fd,fdargs{:});
%end

t=[t;tt+size(p,1)];
p=[p;pp];

bndexpr={'all((p(:,1)-.25).^2+p(:,2).^2<1)','true'};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,[]);
msh.tcurved(:)=true; msh.ecurved(:)=true; msh.ncurved(:)=true;
msh=mshfix(msh);
