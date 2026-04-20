function u=dg2cg(msh,u,W)

if nargin<3
    [r,W]=mapdg2cg(msh.p1);
end

u=W*reshape(permute(u,[1,3,2]),[],size(u,2));
