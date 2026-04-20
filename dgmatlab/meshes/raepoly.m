function prae=raepoly(res)

load raebnd
prae=[x,y];
prae=prae(41:217,:);
prae=interp1(1:size(prae,1),prae,1:res:size(prae,1),'cubic');
