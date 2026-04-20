function perm=mkKperm(msh)

u=msh.p1;
i=reshape(1:prod(size(u)),size(u));
i=permute(i,[1,3,2]);
perm=i(:);
