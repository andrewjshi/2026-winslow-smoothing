function [p,e,t]=squaremesh(m,n,parity)

if nargin<1, m=10; end
if nargin<2, n=m; end
if nargin<3, parity=0; end

% Generate mesh for unit square
[x,y]=ndgrid((0:m-1)/(m-1),(0:n-1)/(n-1));
p=[x(:),y(:)];
if parity==0
  t=int32([1,2,m+2; 1,m+2,m+1]);
else
  t=int32([1,2,m+1; 2,m+2,m+1]);
end
t=kron(t,ones(m-1,1,'int32'))+kron(ones(size(t),'int32'),int32(0:m-2)');
t=kron(t,ones(n-1,1,'int32'))+kron(ones(size(t),'int32'),int32(0:n-2)'*m);
e=int32([1:m,m+1:m:m*n,2*m:m:m*n,m*n-m+2:m*n-1]);

% Reorder triangles in Cartesian order
ix=[];
for i=1:n-1
  ix1=i+(n-1)*(0:m-2);
  ix2=ix1+(n-1)*(m-1);
  if parity==0
    ix12=[ix2;ix1];
  else
    ix12=[ix1;ix2];
  end
  ix=[ix,reshape(ix12,1,[])];
end
t=t(ix,:);
