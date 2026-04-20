function [r,n]=froceval(fro,is,vw)

maken=nargout>=2;

s=fro.vw{is};
nsv=size(s,3);
nsw=size(s,4);

nvw=size(vw,1);
r=zeros(nvw,3);
if maken
  n=zeros(nvw,3);
end

C=[1,0,0,0; 0,0,1,0; -3,3,-2,-1; 2,-2,1,1];
for ivw=1:nvw
  i1v=floor(vw(ivw,1)+1);
  i2v=i1v+1;
  lv=vw(ivw,1)-floor(vw(ivw,1));
  
  i1w=floor(vw(ivw,2)+1);
  i2w=i1w+1;
  lw=vw(ivw,2)-floor(vw(ivw,2));

  if i1v==nsv && lv==0
    i1v=i1v-1;
    i2v=i2v-1;
    lv=1.0;
  end
  if i1v<1 | i2v>nsv, error('Parameter out of range.'); end

  if i1w==nsw && lw==0
    i1w=i1w-1;
    i2w=i2w-1;
    lw=1.0;
  end
  if i1w<1 | i2w>nsw, error('Parameter out of range.'); end
  
  vvec=[1,lv,lv^2,lv^3];
  wvec=[1,lw,lw^2,lw^3];

  xmat=[s(1,1,i1v,i1w), s(1,1,i1v,i2w), s(3,1,i1v,i1w), s(3,1,i1v,i2w);
        s(1,1,i2v,i1w), s(1,1,i2v,i2w), s(3,1,i2v,i1w), s(3,1,i2v,i2w);
        s(2,1,i1v,i1w), s(2,1,i1v,i2w), s(4,1,i1v,i1w), s(4,1,i1v,i2w);
        s(2,1,i2v,i1w), s(2,1,i2v,i2w), s(4,1,i2v,i1w), s(4,1,i2v,i2w)];
  ymat=[s(1,2,i1v,i1w), s(1,2,i1v,i2w), s(3,2,i1v,i1w), s(3,2,i1v,i2w);
        s(1,2,i2v,i1w), s(1,2,i2v,i2w), s(3,2,i2v,i1w), s(3,2,i2v,i2w);
        s(2,2,i1v,i1w), s(2,2,i1v,i2w), s(4,2,i1v,i1w), s(4,2,i1v,i2w);
        s(2,2,i2v,i1w), s(2,2,i2v,i2w), s(4,2,i2v,i1w), s(4,2,i2v,i2w)];
  zmat=[s(1,3,i1v,i1w), s(1,3,i1v,i2w), s(3,3,i1v,i1w), s(3,3,i1v,i2w);
        s(1,3,i2v,i1w), s(1,3,i2v,i2w), s(3,3,i2v,i1w), s(3,3,i2v,i2w);
        s(2,3,i1v,i1w), s(2,3,i1v,i2w), s(4,3,i1v,i1w), s(4,3,i1v,i2w);
        s(2,3,i2v,i1w), s(2,3,i2v,i2w), s(4,3,i2v,i1w), s(4,3,i2v,i2w)];

  r(ivw,1)=vvec*C*xmat*C'*wvec';
  r(ivw,2)=vvec*C*ymat*C'*wvec';
  r(ivw,3)=vvec*C*zmat*C'*wvec';
  
  if maken
    dvvec=[0,1,2*lv,3*lv^2];
    dwvec=[0,1,2*lw,3*lw^2];
    
    dv=[dvvec*C*xmat*C'*wvec';
        dvvec*C*ymat*C'*wvec';
        dvvec*C*zmat*C'*wvec'];
    dw=[vvec*C*xmat*C'*dwvec';
        vvec*C*ymat*C'*dwvec';
        vvec*C*zmat*C'*dwvec'];
    dn=cross(dv,dw);
    dn=dn/norm(dn);
    n(ivw,:)=dn(:)';
  end
end
