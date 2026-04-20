function [x,w]=gaussquad2d(n)

load gaussquad2ddat
nmax=numel(xs);
n=max(n,1);
if n<=numel(xs)
  x=xs{n};
  w=ws{n};
else
  % Map from quad
  [x1,w1]=gaussquad1d(n); x1=2*x1-1; w1=w1*2;
  [x2,y2]=ndgrid(x1,x1);
  x0=[x2(:),y2(:)];
  x=[(1+x0(:,1)-x0(:,2)-x0(:,1).*x0(:,2))/4,(1+x0(:,2))/2];
  w0=w1*w1';
  w=w0(:).*(1-x0(:,2))/4;
end
