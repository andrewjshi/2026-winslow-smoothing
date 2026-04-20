function [x,w]=gaussquad3d(n)

load gaussquad3ddat
nmax=numel(xs);
n=max(n,1);
if n<=numel(xs)
  x=xs{n};
  w=ws{n};
else
  error('Not implemented.');
end
