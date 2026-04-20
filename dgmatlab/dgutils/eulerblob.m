function u0=eulerblob(p,Uinf,p0,r0,d0)

gamma=1.4;
dim=length(Uinf)-2;
pinf=dgeval(Uinf,'p');

blob=exp(-sum((p-repmat(p0(:)',[size(p,1),1,size(p,3)])).^2,2)/r0^2);

rinf=Uinf(1);
uinf=Uinf(2);
u2inf=uinf^2;
if dim>=2
  vinf=Uinf(3);
  u2inf=u2inf+vinf^2;
end
if dim>=3
  winf=Uinf(4);
  u2inf=u2inf+winf^2;
end
r=rinf*(1+d0*blob);
p=pinf*(1+d0*blob);

rE=p/(gamma-1)+1/2*u2inf;
if dim==1
  u0=cat(2,r,r*uinf,rE);
elseif dim==2
  u0=cat(2,r,r*uinf,r*vinf,rE);
elseif dim==3
  u0=cat(2,r,r*uinf,r*vinf,r*winf,rE);
end  
