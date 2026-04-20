function ux=dgdiff(msh,u,comp,chunk)

if nargin<3, comp=1; end
u=u(:,comp,:);
dim=size(msh.p,1);
nt=size(msh.t,2);
ns=size(msh.s,1);

if nargin>=4 & ~isinf(chunk)
    ux=zeros(ns,dim,nt);
    ii=0;
    while ii<nt
        cchunk=min(chunk,nt-ii);
        ix=ii+(1:cchunk);
        msh1=msh;
        msh1.t=msh.t(:,ix);
        msh1.p1=msh.p1(:,:,ix);
        ux(:,:,ix)=dgdiff(msh1,u(:,:,ix));
        ii=ii+cchunk;
    end
    return;
end

switch dim
 case 1
  p1x=permute(msh.p1(:,1,:),[1,3,2]);
  
  gfs=mkvolshapes(double(msh.porder),msh.s,msh.s);
  phiX=permute(gfs(:,2,:),[1,3,2]);
  
  J=phiX'*p1x;
  iJ=1./J;
  
  ns=size(msh.s,1);
  nt=size(msh.t,2);
  ux=zeros(ns,1,nt);
  for it=1:nt
    phix=diag(iJ(:,it))*phiX';
    ux(:,1,it)=phix*u(:,1,it);
  end
 case 2
  p1x=permute(msh.p1(:,1,:),[1,3,2]);
  p1y=permute(msh.p1(:,2,:),[1,3,2]);
  
  gfs=mkvolshapes(double(msh.porder),msh.s,msh.s);
  phiX=permute(gfs(:,2,:),[1,3,2]);
  phiY=permute(gfs(:,3,:),[1,3,2]);
  
  J11=phiX'*p1x;
  J21=phiY'*p1x;
  J12=phiX'*p1y;
  J22=phiY'*p1y;
  J=J11.*J22-J12.*J21;
  
  iJ11=J22./J;
  iJ22=J11./J;
  iJ12=-J12./J;
  iJ21=-J21./J;
  
  ns=size(msh.s,1);
  nt=size(msh.t,2);
  ux=zeros(ns,2,nt);
  for it=1:nt
    phix=diag(iJ11(:,it))*phiX'+diag(iJ12(:,it))*phiY';
    phiy=diag(iJ21(:,it))*phiX'+diag(iJ22(:,it))*phiY';
    
    ux(:,1,it)=phix*u(:,1,it);
    ux(:,2,it)=phiy*u(:,1,it);
  end
 case 3
  p1x=permute(msh.p1(:,1,:),[1,3,2]);
  p1y=permute(msh.p1(:,2,:),[1,3,2]);
  p1z=permute(msh.p1(:,3,:),[1,3,2]);
  
  gfs=mkvolshapes(double(msh.porder),msh.s,msh.s);
  phiX=permute(gfs(:,2,:),[1,3,2]);
  phiY=permute(gfs(:,3,:),[1,3,2]);
  phiZ=permute(gfs(:,4,:),[1,3,2]);
  
  J11=phiX'*p1x;
  J21=phiY'*p1x;
  J31=phiZ'*p1x;
  J12=phiX'*p1y;
  J22=phiY'*p1y;
  J32=phiZ'*p1y;
  J13=phiX'*p1z;
  J23=phiY'*p1z;
  J33=phiZ'*p1z;
  J=J11.*J22.*J33+J12.*J23.*J31+J13.*J21.*J32-J11.*J23.*J32-J12.*J21.*J33-J13.*J22.*J31;

  iJ11=(J22.*J33-J23.*J32)./J;
  iJ12=-(J12.*J33-J13.*J32)./J;
  iJ13=(J12.*J23-J13.*J22)./J;
  iJ21=-(J21.*J33-J23.*J31)./J;
  iJ22=(J11.*J33-J13.*J31)./J;
  iJ23=-(J11.*J23-J13.*J21)./J;
  iJ31=(J21.*J32-J22.*J31)./J;
  iJ32=-(J11.*J32-J12.*J31)./J;
  iJ33=(J11.*J22-J12.*J21)./J;
  
  ns=size(msh.s,1);
  nt=size(msh.t,2);
  ux=zeros(ns,3,nt);
  for it=1:nt
    phix=diag(iJ11(:,it))*phiX'+diag(iJ12(:,it))*phiY'+diag(iJ13(:,it))*phiZ';
    phiy=diag(iJ21(:,it))*phiX'+diag(iJ22(:,it))*phiY'+diag(iJ23(:,it))*phiZ';
    phiz=diag(iJ31(:,it))*phiX'+diag(iJ32(:,it))*phiY'+diag(iJ33(:,it))*phiZ';
    
    ux(:,1,it)=phix*u(:,1,it);
    ux(:,2,it)=phiy*u(:,1,it);
    ux(:,3,it)=phiz*u(:,1,it);
  end
end
