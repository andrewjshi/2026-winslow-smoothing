function [KWP,KWPsizes,KWPegix]=mkKWprojections(dim,porder)

ns=nchoosek(porder+dim,dim);
KWP=zeros(ns,ns,porder,porder);
KWPsizes=ones(porder+1,1,'int32');
for pH=0:porder
  KWPsizes(pH+1,1)=nchoosek(pH+dim,dim);
  KWPsizes(pH+1,2)=nchoosek(pH+dim-1,dim-1);
end

for pH=1:porder
  for pL=0:pH-1
    sH=lagrangepnts(pH,dim);
    VH=kwvander(pH,2*sH(:,1:dim)-1);
    sL=lagrangepnts(pL,dim);
    VL=kwvander(pL,2*sL(:,1:dim)-1);
    
    VH=VH(:,1:size(VL,1));
    P=VH/VL;
    
    KWP(1:size(P,1),1:size(P,2),pH,pL+1)=P;
  end
end

KWPegix=zeros(porder+1,dim+1,porder+1,'int32');
for p=0:porder
  egix=mkegix(p,dim);
  KWPegix(1:size(egix,1),:,p+1)=egix;
end

function egix=mkegix(porder,dim)

s=lagrangepnts(porder,dim);
sbnd=lagrangepnts(porder,dim-1);

egix=zeros(size(sbnd,1),0,'int32');
facemap=mkfacemap(dim);
for ii=1:dim+1
  ic=facemap(:,ii);
  [foo,ix]=ismember(round(porder*sbnd),round(porder*s(:,ic)),'rows');
  egix=[egix,int32(ix-1)];
end
