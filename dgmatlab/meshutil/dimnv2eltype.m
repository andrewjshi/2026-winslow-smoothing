function eltype=dimnv2eltype(dim,nv)

if nv==dim+1
    eltype=t_simplex;
elseif nv==2^dim
    eltype=t_block;
else
    error('Unknown element type');
end
