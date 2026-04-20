function u=cg2dg(msh,u,W1)

if nargin<3
    [r,W,W1]=mapdg2cg(msh.p1);
end

u=permute(reshape(W1'*u, size(msh.p1,1),size(msh.p1,3),[]),[1,3,2]);
