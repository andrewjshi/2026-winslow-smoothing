function msh=mshdiamond2(nref,hs,box)
%MSHDIAMOND2 Unstructured mesh for diamond geometry.
%   msh=mshdiamond2(nref,hmax)

if nargin<1, nref=4; end
if nargin<2, hs=[.2,2]; end
if nargin<3, box=[-2.5,7,-2.5,2.5]; end

pv1=[1,0; 0,0.5; -1,0; 0,-.5; 1,0];
for iref=1:nref
  mid=(pv1(1:end-1,:)+pv1(2:end,:))/2;
  pvnew=zeros(size(pv1,1)+size(mid,1),2);
  pvnew(1:2:end,:)=pv1;
  pvnew(2:2:end,:)=mid;
  pv1=pvnew;
end

pv2=[box(1),box(3); box(2),box(3); box(2),box(4); box(1),box(4)];

[p,t]=polymesh({pv1,pv2},[0,1],[0,1;1,0],[hs(2),1.6],@href,hs);

% Create boundary numbering and msh structure
bndexpr={'all(abs(p(:,1))<1.1 & abs(p(:,2))<1.1)', ...
         'any(abs(p(:,1))>1.1 | abs(p(:,2))>1.1)'};
msh=ml2msh(p,t,bndexpr);
msh=mshcurved(msh,[]);

function h=href(p,hs)

pv1=[0,-.75; 7,-1.5; 7,1.5; 0,.75];
inside1=inpolygon(p(:,1),p(:,2),pv1(:,1),pv1(:,2));
h=inf*ones(size(p,1),1);
h(inside1)=hs(1);
