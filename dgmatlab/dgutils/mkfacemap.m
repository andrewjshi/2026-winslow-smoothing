function facemap=mkfacemap(dim,eltype)
%MKFACEMAP  Indices of faces in element

if nargin<2, eltype=t_simplex; end

% Should match corresponding definition in src/dgutil.cpp

switch eltype
  case t_simplex
    switch dim
      case 1
        facemap=[2,1];
      case 2
        facemap=[2,3;3,1;1,2]';
      case 3
        facemap=[2,3,4;1,4,3;4,1,2;3,2,1]';
      otherwise
        error('Dimension not implemented.');
    end
  case t_block
    switch dim
      case 1
        facemap=[1,2];
      case 2
        facemap=[3,1;2,4;1,2;4,3]';
      case 3
        facemap=[3,1,7,5; 2,4,6,8; 1,2,5,6; 4,3,8,7; 3,4,1,2; 5,6,7,8]';
      otherwise
        error('Dimension not implemented.');
    end
  otherwise
    error('Element type not implemented.');
end
