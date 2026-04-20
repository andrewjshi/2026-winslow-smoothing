function msh=mshfoil1(R,alpha,nvals)

if nargin<1, R=10; end
if nargin<2, alpha=3; end
if nargin<3, nvals=[1,3,100,2,10,0.5]; end

load foil1bnd
p1=p;
np1=size(p1,1);

pp1=interp1(1:np1,p1,1:0.25:np1,'cubic');

fd=inline('ddiff(dcircle(p,0.25,0,R),dpoly(p,pp))','p','pp','R');
fdargs={pp1,R};

cn=cos(alpha*pi/180);
sn=sin(alpha*pi/180);
l=-.75*cn+sqrt(.75^2*cn^2-.75^2+R^2);
s=block(nvals(3),nvals(4),[0;l*nvals(6)])';
pwake=[1+cn*s, sn*s];

beta=atan2(pwake(end,2),pwake(end,1)-.25);

n=16;
phi=beta+2*pi*(1:n)'/n;
pc=[.25+R*cos(phi),R*sin(phi)];

s=0:np1-1;
d=round(nvals(1)+(s-(np1-1)/2).^2*nvals(2)/((np1-1)/2)^2);

i=1;
ix=1;
while ix(i)<np1-1
  ix(i+1)=ix(i)+d(ix(i));
  i=i+1;
end
ix(end)=[];

[p,t]=polymesh({p1(ix,:),pc,pwake},[1,1,0],[0,1;1,0;1,1],[nvals(5),1.8]);
for ii=1:5
  p=bndproj(p,t,fd,fdargs{:});
end

bndexpr={'all((p(:,1)-.25).^2+p(:,2).^2<1)','true'};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,[1,2]);
msh=mshfix(msh);
