function dim=nv2dim(nv,eltype)

switch eltype
  case t_simplex
    dim=nv-1;
  case t_block
    dim=log2(nv);
  otherwise
    error('Element type not implemented.');
end
