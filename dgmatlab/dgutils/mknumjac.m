function J=mknumjac(u,msh,data,phys,resname)

R=feval(resname,u,msh,data,phys);
[ns,nc,nt]=size(R);
N=ns*nc*nt;
nu=size(u,2);
J=sparse(N,N);
for ii=1:ns
  for jj=1:nc
    for kk=1:nt
      u1=u;
      u1(ii,jj,kk)=u(ii,jj,kk)+1e-6;
      R1=feval(resname,u1,msh,data,phys);
      u1(ii,jj,kk)=u(ii,jj,kk)-1e-6;
      R2=feval(resname,u1,msh,data,phys);
      J(:,ii+(jj-1)*ns+(kk-1)*ns*nc)=(R1(:)-R2(:))/2e-6;
    end
  end
end
