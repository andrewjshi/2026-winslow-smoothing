function msh=mshbndreorder(msh,perm)

if iscell(perm) % New syntax
    for i=1:numel(perm)
        ix = find(ismember(msh.t2t, -perm{i}));
        msh.t2t(ix) = -i;
    end
    msh.bndexpr = {};
else
    [foo,invperm]=sort(perm);
    ix=find(msh.t2t<0);
    msh.t2t(ix)=-invperm(-msh.t2t(ix));
    msh.bndexpr=msh.bndexpr(perm);
end
