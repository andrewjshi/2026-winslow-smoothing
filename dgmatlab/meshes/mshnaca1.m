function msh=mshnaca1(angle,hs)

if nargin<1, angle=5; end
if nargin<2, hs=[0.025,0.1,0.1,0.2,0.4,1]; end

xlim=[0,1.0089304129];
a=.12/.2*[0.2969,-0.1260,-0.3516,0.2843,-0.1015];
theta=angle*pi/180;

fd=@dnaca;
fdargs={a,theta};

x1=(hs(1)/a(1))^2;
x2=xlim(2)-hs(2);
x=[xlim(1),x1,x2,xlim(2)]';

[x,y]=limit_arcratio(x,hs(3),a);
R=[cos(theta),-sin(theta); sin(theta),cos(theta)];

pv1=[x(end:-1:1),y(end:-1:1); x(2:end-1),-y(2:end-1)]*R;
pv2=[-1,-2; 4,-2; 4,2; -1,2];

[p,t]=polymesh({pv1,pv2},[1,1],[0,1;1,0],[hs(end),1.5],@href,hs);

fdbox=@(p) ddiff(drectangle(p,-1,4,-2,2),dnaca(p,a,theta));
for ii=1:5
  p=bndproj(p,t,fdbox);
end

bndexpr={'all(p(:,2)<-2+1e-6)','all(p(:,1)>4-1e-6)', ...
         'all(p(:,2)>2-1e-6)','all(p(:,1)<-1+1e-6)', ...
         'all((p(:,1)-.5).^2+p(:,2).^2<1)'};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,5);
msh=mshfix(msh);

function h=href(p,hs)

pv1=[4,.75;0,.25;0,-.25;4,-.75;4,.75];
pv2=[2,.5;0,.25;0,-.25;2,-.5;2,.5];
inside1=inpolygon(p(:,1),p(:,2),pv1(:,1),pv1(:,2));
inside2=inpolygon(p(:,1),p(:,2),pv2(:,1),pv2(:,2));
h=inf*ones(size(p,1),1);
h(inside1)=hs(5);
h(inside2)=hs(4);

function [x,y]=limit_arcratio(x,hmax,a)

while 1
  y=polyval([a(5:-1:2),0],x)+a(1)*sqrt(x);
  s=sqrt((x(2:end)-x(1:end-1)).^2+(y(2:end)-y(1:end-1)).^2);
  
  [maxs,ix]=max(s);
  if maxs>1.4*hmax
    x=[x(1:ix);(x(ix)+x(ix+1))/2;x(ix+1:end)];
    continue;
  end
  
  ratio=s(2:end)./s(1:end-1);
  ix=find(ratio>2 | ratio<1/2);
  
  if isempty(ix)
    break;
  else
    ix=ix(1);
    if ratio(ix)>2
      x=[x(1:ix+1);(x(ix+1)+x(ix+2))/2;x(ix+2:end)];
    else
      x=[x(1:ix);(x(ix)+x(ix+1))/2;x(ix+1:end)];
    end
  end
end
