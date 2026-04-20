function [blkP,blkPsizes,blkPegix]=mkblockprojections(dim,porder)

ns=(porder+1)^dim;
blkP=zeros(ns,ns,porder,porder);
blkPsizes=ones(porder+1,2,'int32');
for pH=0:porder
  blkPsizes(pH+1,1)=(pH+1)^dim;
  blkPsizes(pH+1,2)=(pH+1)^(dim-1);
end

for pH=1:porder
  for pL=0:pH-1
    V0H=plegendre(lobattopoints(pH+1),pH);
    nH=size(V0H,1);
    V0L=plegendre(lobattopoints(pL+1),pL);
    nL=size(V0L,1);
    P0=V0H(:,1:nL)/V0L;

    switch dim
      case 1
        P=P0;
      case 2
        P=reshape(permute(reshape(P0(:)*P0(:)',[nH,nL,nH,nL]),[1,3,2,4]),[nH^2,nL^2]);
      case 3
        P=reshape(permute(reshape(P0(:)*P0(:)',[nH,nL,nH,nL]),[1,3,2,4]),[nH^2,nL^2]);
        P=reshape(permute(reshape(P(:)*P0(:)',[nH^2,nL^2,nH,nL]),[1,3,2,4]),[nH^3,nL^3]);
      otherwise
        error('Dimension not implemented.');
    end
    blkP(1:size(P,1),1:size(P,2),pH,pL+1)=P;
  end
end

blkPegix=zeros(porder+1,2*dim,porder+1,'int32');
for p=0:porder
  egix=mkegix(p,dim);
  blkPegix(1:size(egix,1),:,p+1)=egix;
end

function egix=mkegix(porder,dim)

if porder==0
    egix=zeros(1,2*dim,'int32');
    return
end

[s,foo,sbnd,foo]=lagrangepnts(porder,dim,t_block);
snap=@(x) round(x*1e8)/1e8;

egix=zeros(size(sbnd,1),0,'int32');
facemap=mkfacemap(dim,t_block);
for ii=1:size(facemap,2)
    ic=facemap(:,ii);
    [foo,ix]=ismember(snap(sbnd),snap(s(:,ic)),'rows');
    egix=[egix,int32(ix-1)];
end
