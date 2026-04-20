function u=freestream(msh,phys)

u=repmat(phys.uinf,[size(msh.s,1),1,size(msh.t,2)]);
