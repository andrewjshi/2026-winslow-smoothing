function msh=fro2msh(fro,p,t)
%FRO2MSH
%
%   MSH=FRO2MSH(FRO,P,T)
%
% Given surface mesh and definition FRO + tetrahedral
% mesh P,T, create 3DG msh structure MSH.

[t2t,t2n]=mkt2t(t);

fromap=0*t;
facemap=mkfacemap(3);
for j=1:4
  faces=t(:,facemap(:,j));
  [bnd,loc]=ismember(sort(faces,2),sort(fro.e(:,1:3),2),'rows');
  t2t(bnd,j)=-fro.e(loc(bnd),4)+1;
  fromap(bnd,j)=loc(bnd);
end

msh.p=p';
msh.t=int32(t'-1);
msh.t2t=int32(t2t'-1);
msh.t2n=int32(t2n'-1);
msh=mshcurved(msh,[]);
msh.eltype = int32(0);
msh.fromap=fromap;
msh.fro=fro;
