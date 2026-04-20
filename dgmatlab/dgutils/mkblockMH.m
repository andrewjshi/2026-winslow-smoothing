function MH=mkblockMH(M,s0,dim)

ns=numel(s0);
porder=ns-1;

V0=plegendre(2*s0-1,porder);
PL=zeros(ns,ns);
PL(1:ns-1,1:ns-1)=eye(ns-1);
FL0=V0*PL*inv(V0);

switch dim
  case 1
    FL=FL0;
  case 2
    FL=reshape(permute(reshape(FL0(:)*FL0(:)',[ns,ns,ns,ns]),[1,3,2,4]),[ns^2,ns^2]);
  case 3
    FL=reshape(permute(reshape(FL0(:)*FL0(:)',[ns,ns,ns,ns]),[1,3,2,4]),[ns^2,ns^2]);
    FL=reshape(permute(reshape(FL(:)*FL0(:)',[ns^2,ns^2,ns,ns]),[1,3,2,4]),[ns^3,ns^3]);
  otherwise
    error('Dimension not implemented.');
end

FH=eye(size(FL))-FL;
MH=FH'*M*FH;
