function s=sdistr(NX,NY,L,ws,mins)

if length(L)==1, L=L*ones(1,ws(1)); end
L=L(1:ws(1));

ds0=logdistr(mins(1),L(end),NY-1);
ds=repmat(ds0,NX,1);

for ii=1:ws(1)-1
  wt=ii/(ws(1)-1);
  cmins=exp(log(mins(2))*wt+log(mins(3))*(1-wt));
  cds=ds0;
  ix=max(find(cds<cmins));
  cds(1:ix)=cmins;
  cds=cds/sum(cds)*L(ii);
  ds(ii,:)=cds;
  ds(NX-ii+1,:)=cds;
end
s=[zeros(NX,1),cumsum(ds,2)];

function ds=logdistr(mins,L,NY)

pol=zeros(NY+1,1);
pol(1)=mins;
pol(end-1)=-L;
pol(end)=-(mins-L);
r=roots(pol);
r=max(r(imag(r)==0));

ix=0:(NY-1);
ds=mins*r.^ix;
