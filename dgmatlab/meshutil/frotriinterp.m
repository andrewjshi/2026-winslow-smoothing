function [r,n]=frotriinterp(fro,ie,s)

maken=nargout>=2;

s=s(:,1:2);
s=[s,1-sum(s,2)];
ce=fro.e(ie,:);
ns=size(s,1);

r=zeros(ns,3);
if maken
  n=zeros(ns,3);
end
for is=1:ns
  cs=s(is,:);

  triinterp=true;
  if sum(cs==0)==2 % Corner
    r(is,:)=cs*fro.p(ce(1:3),:);
    triinterp=false;
  elseif sum(cs==0)==1 % Edge
    ix=find(cs~=0);
    r(is,:)=froedgeinterp(fro,ce(ix),cs(ix),true);
    if ~any(isnan(r(is,:))) % On a curve
      triinterp=false;
    end
  end

  if triinterp || maken % General case, interpolate using surface
    ips=ce(4);
    cps=fro.ps{ips};
    
    vws=zeros(3,2);
    vws(1,:)=cps(cps(:,1)==ce(1),2:3);
    vws(2,:)=cps(cps(:,1)==ce(2),2:3);
    vws(3,:)=cps(cps(:,1)==ce(3),2:3);
    
    vw=cs*vws;

    if maken
      [cr,cn]=froseval(fro,ips,vw);
      n(is,:)=cn;
    else
      cr=froseval(fro,ips,vw);
    end
    if triinterp
      r(is,:)=cr;
    end
  end
end
