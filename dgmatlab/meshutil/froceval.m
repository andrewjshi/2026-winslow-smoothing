function [r,t]=froceval(fro,ic,u)

maket=nargout>=2;

c=fro.u{ic};
nc=size(c,3);

nu=numel(u);
r=zeros(nu,3);
if maket
  t=zeros(nu,3);
end

C=[1,0,0,0; 0,0,1,0; -3,3,-2,-1; 2,-2,1,1];
for iu=1:nu
  i1=floor(u(iu)+1);
  i2=i1+1;
  lu=u(iu)-floor(u(iu));
  
  if i1==nc && lu==0
    i1=i1-1;
    i2=i2-1;
    lu=1.0;
  end
  
  if i1<1 | i2>nc, error('Parameter out of range.'); end

  uvec=[1,lu,lu^2,lu^3];
  rmat=[c(1,:,i1); c(1,:,i2); c(2,:,i1); c(2,:,i2)];
  
  r(iu,:)=uvec*C*rmat;
  
  if maket
    duvec=[0,1,2*lu,3*lu^2];
    t(iu,:)=duvec*C*rmat;
    t(iu,:)=t(iu,:)/norm(t(iu,:));
  end
end
