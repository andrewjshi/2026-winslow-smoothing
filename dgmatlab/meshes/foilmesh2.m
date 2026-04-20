function [p,t]=foilmesh2(res,nbndlayer,ratiobndlayer,onlybndlayer,doplot)

% res = [x1,x2,x3,y1,ybnd]

if nargin<1 | isempty(res), res=[4,4,4,2,2]; end
if nargin<2 | isempty(nbndlayer), nbndlayer=11; end
if nargin<3 | isempty(ratiobndlayer), ratiobndlayer=5; end
if nargin<4 | isempty(onlybndlayer), onlybndlayer=0; end
if nargin<5 | isempty(doplot), doplot=0; end

load grid.ht13
x=grid(:,1);
y=grid(:,2);
ix=find(diff(x)<0);
n=ix(1);
x=reshape(x,n,[]);
y=reshape(y,n,[]);

load bl.ht13
xb=bl(:,1);
yb=bl(:,2);

bix=find(diff(xb)<0);
bix=bix(end);

xb1=xb(1:bix);
yb1=yb(1:bix);
xb2=xb(bix+1:2*bix);
yb2=yb(bix+1:2*bix);
xb3=xb(2*bix+1:end);
yb3=yb(2*bix+1:end);

dy=diff(y,[],2);
[foo,ix]=min(dy(1,:));
i1=find(dy(:,ix)>1e-8);
i1=i1(1)-1;

x1=x(1:i1,[1:ix,ix+2:end]);
y1=y(1:i1,[1:ix,ix+2:end]);
x2=x(i1:i1+bix-1,1:ix);
y2=y(i1:i1+bix-1,1:ix);
x3=x(i1:i1+bix-1,ix+1:end);
y3=y(i1:i1+bix-1,ix+1:end);
x4=x(i1+bix-1:end,1:ix);
y4=y(i1+bix-1:end,1:ix);
x5=x(i1+bix-1:end,ix+1:end);
y5=y(i1+bix-1:end,ix+1:end);

xb3=[xb1(end); xb3; (x4(end,end)+x5(end,1))/2];
yb3=[yb1(end); yb3; (y4(end,end)+y5(end,1))/2];

s=block(nbndlayer,ratiobndlayer,[0;1]);
x6=xb2*(1-s)+x2(:,end)*s;
y6=yb2*(1-s)+y2(:,end)*s;
x7=xb1*(1-s)+x3(:,1)*s;
y7=yb1*(1-s)+y3(:,1)*s;
x8=xb3*(1-s)+x4(:,end)*s;
y8=yb3*(1-s)+y4(:,end)*s;
x9=xb3*(1-s)+x5(:,1)*s;
y9=yb3*(1-s)+y5(:,1)*s;

xblks={x1,x2,x3,x4,x5,x6,x7,x8,x9};
yblks={y1,y2,y3,y4,y5,y6,y7,y8,y9};
for i=1:9
  blks{i}=cat(3,xblks{i},yblks{i});
end

dd=[res(1),res(4);
    res(2),res(4);
    res(2),res(4);
    res(3),res(4);
    res(3),res(4);
    res(2),res(5);
    res(2),res(5);
    res(3),res(5);
    res(3),res(5)];

for i=1:length(blks)
  blks{i}=blks{i}(1:dd(i,1):end,1:dd(i,2):end,:);
end

blk0=1;
if onlybndlayer
  blk0=6;
end

if doplot
  col='rbmcygrrg';
  clf,axis equal,view(2)
  for i=blk0:length(blks)
    surface(blks{i}(:,:,1),blks{i}(:,:,2),0*blks{i}(:,:,1), ...
            'edgecolor','k','facecolor',col(i));
  end
end

p=zeros(0,2);
t=zeros(0,3);
for i=blk0:length(blks)
  blk=blks{i};
  pp=reshape(blk,[],2);
  [foo,foo,tt]=squaremesh(size(blk,1),size(blk,2),1);
%  [pp,tt]=block2pt(permute(blks{i},[3,1,2]));
  t=[t;tt+size(p,1)];
  p=[p;pp];
end
[p,t]=fixmesh(p,t);
