function MH=mkMH(porder,M,s)

ns=size(s,1);
dim=size(s,2)-1;

maxloworder=porder-1;
if maxloworder>=0
  ns1=nchoosek(maxloworder+dim,dim);
else
  ns1=1;
end

V=kwvander(porder,2*s(:,1:dim)-1);
PH=diag([zeros(1,ns1),ones(1,ns-ns1)]);
FH=V*PH*inv(V);
MH=FH'*M*FH;
