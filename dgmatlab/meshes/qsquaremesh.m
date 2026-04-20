function [p,e,q]=qsquaremesh(m,n)

if nargin<1, m=10; end
if nargin<2, n=m; end

% Generate mesh for unit square
[x,y]=ndgrid((0:m-1)/(m-1),(0:n-1)/(n-1));
p=[x(:),y(:)];
q=int32([1,2,m+1,m+2]);
q=kron(q,ones(m-1,1,'int32'))+kron(ones(size(q),'int32'),int32(0:m-2)');
q=kron(q,ones(n-1,1,'int32'))+kron(ones(size(q),'int32'),int32(0:n-2)'*m);
e=int32([1:m,m+1:m:m*n,2*m:m:m*n,m*n-m+2:m*n-1]);
