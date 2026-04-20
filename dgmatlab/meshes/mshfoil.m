function [msh,pdist]=mshfoil(foilp,R,alpha,nvals)

if nargin<2, R=10; end
if nargin<3, alpha=0; end
if nargin<4, nvals=[50,5,20,10,5]; end

% nvals = [nfoil, ratio_foil, nwake, ratio_wake, hmax]

if isnumeric(foilp)
elseif ischar(foilp)
  switch foilp
   case 'ht13'
    load foil1bnd
    foilp=p;
  end
end

p1=foilp;
np1=size(p1,1);

cn=cos(alpha*pi/180);
sn=sin(alpha*pi/180);
l=-.75*cn+sqrt(.75^2*cn^2-.75^2+R^2); % This might assume chord = 1 ...
s=block(nvals(3),nvals(4),[0;l*1.0])';
pwake=[max(p(:,1))+cn*s, sn*s];

beta=atan2(pwake(end,2),pwake(end,1)-.25);

n=16;
phi=beta+2*pi*(1:n)'/n;
pc=[.25+R*cos(phi),R*sin(phi)];

arcs=[0;cumsum(sqrt(sum(diff(p1,[],1).^2,2)))];
ss=block(nvals(1),-1/nvals(2),[0;arcs(end)]);
p0=interp1(arcs,p1,ss,'cubic');
p0=p0(1:end-1,:);

[p,t]=polymesh({p0,pc,pwake},[1,1,0],[0,1;1,0;1,1],[nvals(5),1.8]);

pp1=interp1(1:np1,p1,1:0.25:np1,'cubic');
fd=inline('ddiff(dcircle(p,0.25,0,R),dpoly(p,pp))','p','pp','R');
fdargs={pp1,R};
for ii=1:5
  p=bndproj(p,t,fd,fdargs{:});
end

bndexpr={'all((p(:,1)-.25).^2+p(:,2).^2<1)','true'};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,[1,2]);
msh=mshfix(msh);

pdist=[pp1;pwake];
