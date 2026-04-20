function [x,w]=gaussquad(porder,dim)

if nargin<2, dim=1; end

switch dim
 case 0
  x=zeros(1,0); w=1;
 case 1
  [x,w]=gaussquad1d(porder);
 case 2
  [x,w]=gaussquad2d(porder);
 case 3
  [x,w]=gaussquad3d(porder);
 otherwise
  error('Dimension not implemented.');
end
