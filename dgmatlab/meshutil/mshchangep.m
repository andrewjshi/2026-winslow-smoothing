function msh1 = mshchangep(msh, porder)
%MSHCHANGEP  Change polynomial degree of high-order 3DG msh by interpolation

msh1 = nodealloc(msh, porder);
if msh.eltype == t_simplex
    msh1.p1 = dginterp(msh.p1, msh.porder, msh1.porder, size(msh.p,1), msh.eltype);
elseif msh.eltype == t_block
    msh1.p1 = dginterp(msh.p1, msh.s0, msh1.s0, size(msh.p,1), msh.eltype);
end

