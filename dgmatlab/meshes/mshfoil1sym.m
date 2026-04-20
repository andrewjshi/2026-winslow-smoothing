function msh=mshfoil1sym(R,nvals,nref,hs)

if nargin<1, R=10; end
if nargin<2, nvals=[1,3,20,2,10,0.5,10,5]; end
if nargin<3, nref=0; end
if nargin<4, hs=[2,2]; end

load foil1bnd
p1orig=p;
p1=p(1:179,:);
p1(end,:)=[0,0];
np1=size(p1,1);

n=8;
phi=pi*(0:n)'/n;
pc=[R*cos(phi),R*sin(phi)];

s=block(nvals(3),nvals(4),nvals(6)*[0;R-1])';
pwake=[1+s, 0*s];

pp1=interp1(1:np1,p1,1:0.25:np1,'cubic');

fd=inline('ddiff(ddiff(dcircle(p,0,0,R),p(:,2)),dpoly(p,pp))','p','pp','R');
fdargs={pp1,R};

p1=flipud(p1);
s=block(nvals(7),nvals(8),[0;1]);

p1=interp1(1:np1,p1,1+s*(np1-1),'cubic');
pv=[p1;pwake(2:end,:);pc];

[p,t]=polymesh({pv},[1],[1,0],[nvals(5),1.8],@href,hs);
for ii=1:5
  p=bndproj(p,t,fd,fdargs{:});
end

for iref=1:nref
  [p,t]=uniref(p,t);
  for ii=1:5
    p=bndproj(p,t,fd,fdargs{:});
  end
end

t=[t;t+size(p,1)];
p=[p;p(:,1),-p(:,2)];
[p,t]=fixmesh(p,t,1e-6);

fd=inline('ddiff(dcircle(p,0,0,R),dpoly(p,pp))','p','pp','R');
fdargs={p1orig,R};

bndexpr={'all((p(:,1)-.25).^2+p(:,2).^2<1)','true'};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,[1,2]);
msh=mshfix(msh);

function h=href(p,hs)

pv1=[-.1,0; 1,-1; 1,1; -.1,0];
pv2=[1,-1;3,-1;3,1;1,1;1,-1];
inside1=inpolygon(p(:,1),p(:,2),pv1(:,1),pv1(:,2));
inside2=inpolygon(p(:,1),p(:,2),pv2(:,1),pv2(:,2));
h=inf*ones(size(p,1),1);
h(inside1)=hs(1);
h(inside2)=hs(2);
