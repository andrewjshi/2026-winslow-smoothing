function data=dginit(msh,pars)
%DGINIT  Initialize DG data structures (any dimension)
%
%    Syntax: data=dginit(msh)
%
%    msh:   mesh structure
%    data:  data structure

if nargin<2, pars=struct; end
if ~isfield(pars,'gaussfac'), pars.gaussfac=3; end
if ~isfield(pars,'gaussadds'), pars.gaussadds=false; end

porder=double(msh.porder);
s=msh.s;
sbnd=msh.sbnd;
nsbnd=size(msh.sbnd,1);
dim=size(msh.p,1);

switch msh.eltype
  case t_simplex
    [gx,gw]=gaussquad(pars.gaussfac*porder,dim);
    gx=[gx,1-sum(gx,2)];
    if pars.gaussadds, gx=[gx;s]; gw=[gw;0*s(:,1)]; end
    gfs=mkvolshapes(porder,s,gx);
    wgfs=gfs.*repmat(permute(gw,[2,3,1]),[size(gfs,1),dim+1]);
    [lMrefR,lMref]=mkcholM0(porder,s);
    lMrefH=mkMH(porder,lMref,s);

    [egx,egw]=gaussquad(pars.gaussfac*porder,dim-1);
    egx=[egx,1-sum(egx,2)];
    if pars.gaussadds, egx=[egx;sbnd]; egw=[egw;0*sbnd(:,1)]; end
    efs=mkvolshapes(porder,sbnd,egx);
    wefs=efs.*repmat(permute(egw,[2,3,1]),[size(efs,1),dim]);
    gx0 = []; gw0 = []; gfs0 = []; wgfs0 = []; Nx0 = [];
    lM0 = []; lM0R = []; lM0H = [];
    
    [KWP,KWPsizes,KWPegix]=mkKWprojections(dim,porder);
  case t_block
    gx0=gaussquad(max(pars.gaussfac*porder,1));
    [gx,gw,gfs,wgfs]=mkblockshapefcns(dim,msh.s0,gx0);
    [egx,egw,efs,wefs]=mkblockshapefcns(dim-1,msh.s0,gx0);
    [gx0,gw0,gfs0,wgfs0]=mkblockshapefcns(1,msh.s0,gx0);
    Nx0 = mkNx0(msh.s0);
    
    lMref=permute(gfs(:,1,:),[1,3,2])*diag(gw,0)* ...
        permute(gfs(:,1,:),[1,3,2])';
    lMrefR=chol(lMref);
    lMrefH=mkblockMH(lMref,msh.s0,dim);

    lM0=permute(gfs0(:,1,:),[1,3,2])*diag(gw0,0)* ...
        permute(gfs0(:,1,:),[1,3,2])';
    lM0R=chol(lM0);
    lM0H=mkblockMH(lM0,msh.s0,1);

    [KWP,KWPsizes,KWPegix]=mkblockprojections(dim,porder);
  otherwise
    error('Unknown element type');
end

snap=@(x) round(x*1e8)/1e8;

egix=zeros(nsbnd,0,'int32');
facemap=mkfacemap(dim,msh.eltype);
for ii=1:size(facemap,2)
    ic=facemap(:,ii);
    %    Used to assume equi-distant points
    %    [foo,ix]=ismember(round(porder*sbnd),round(porder*s(:,ic)),'rows');
    [foo,ix]=ismember(snap(porder*sbnd),snap(porder*s(:,ic)),'rows');
    egix=[egix,int32(ix-1)];
end

permmap=zeros(nsbnd,0,'int32');
perms={1,[2,1],[1,3,2;2,1,3;3,2,1];
       1,[2,1],[1,3,2,4; 2,1,4,3; 3,4,1,2; 4,2,3,1]};

permegix=zeros(nsbnd,size(facemap,2),0,'int32');
for ii=1:size(perms{msh.eltype+1,dim},1)
    perm=perms{msh.eltype+1,dim}(ii,:);
    %    Used to assume equi-distant points
    %    [foo,ix]=ismember(round(porder*sbnd),round(porder*sbnd(:,perm)),'rows');
    [foo,ix]=ismember(snap(porder*sbnd),snap(porder*sbnd(:,perm)),'rows');
    permmap=[permmap,int32(ix-1)];
    permegix=cat(3,permegix,egix(ix,:));
end

data=struct('eltype',int32(msh.eltype), ...
            'gx0',gx0,'gw0',gw0,'gfs0',gfs0,'wgfs0',wgfs0, 'Nx0', Nx0, ...
            'gx',gx,'gw',gw,'gfs',gfs,'wgfs',wgfs,'lMref',lMref,'lMrefR',lMrefR,'lMrefH',lMrefH, ...
            'lM0',lM0','lM0R',lM0R,'lM0H',lM0H, ...
            'egx',egx,'egw',egw,'efs',efs,'wefs',wefs, ...
            'egix',egix,'permmap',permmap,'permegix',permegix, ...
            'KWP',KWP,'KWPsizes',KWPsizes,'KWPegix',KWPegix);

function Nx0 = mkNx0(s0)

ns=numel(s0);
porder=ns-1;

[PHI,PHIx]=plegendre(2*s0-1,porder);
Nx0 = 2 * (PHIx / PHI);

function [gx,gw,gfs,wgfs]=mkblockshapefcns(dim,s0,gx0)

ns=numel(s0);
porder=ns-1;
ng=size(gx0,1);

[PHI,PHIx]=plegendre(2*s0-1,porder);
PHIx=PHIx*2;
[gN,gNx]=plegendre(2*gx0-1,porder);
gNx=gNx*2;
gN=gN/PHI;
gNx=gNx/PHI;

gw0=gaussweights(2*gx0-1)/2;
gwN=diag(gw0)*gN;
gwNx=diag(gw0)*gNx;

switch dim
  case 0
    gx=1;
    gw=1;
    gfs=1;
  case 1
    gx=gx0;
    gw=gw0;
    gfs=permute(cat(3,gN,gNx),[2,3,1]);
  case 2
    [gx1,gx2]=ndgrid(gx0,gx0);
    gx=[gx1(:),gx2(:)];
    gw=reshape(gw0(:)*gw0(:)',[],1);
    gNN=reshape(permute(reshape(gN(:)*gN(:)',[ng,ns,ng,ns]),[1,3,2,4]),[ng^2,ns^2]);
    gNNx=reshape(permute(reshape(gNx(:)*gN(:)',[ng,ns,ng,ns]),[1,3,2,4]),[ng^2,ns^2]);
    gNNy=reshape(permute(reshape(gN(:)*gNx(:)',[ng,ns,ng,ns]),[1,3,2,4]),[ng^2,ns^2]);
    gfs=permute(cat(3,gNN,gNNx,gNNy),[2,3,1]);
  case 3
    [gx1,gx2,gx3]=ndgrid(gx0,gx0,gx0);
    gx=[gx1(:),gx2(:),gx3(:)];
    gw=reshape(reshape(gw0(:)*gw0(:)',[],1)*gw0(:)',[],1);
    gNN=reshape(permute(reshape(gN(:)*gN(:)',[ng,ns,ng,ns]),[1,3,2,4]),[ng^2,ns^2]);
    gNNx=reshape(permute(reshape(gNx(:)*gN(:)',[ng,ns,ng,ns]),[1,3,2,4]),[ng^2,ns^2]);
    gNNy=reshape(permute(reshape(gN(:)*gNx(:)',[ng,ns,ng,ns]),[1,3,2,4]),[ng^2,ns^2]);
    gNNN=reshape(permute(reshape(gNN(:)*gN(:)',[ng^2,ns^2,ng,ns]),[1,3,2,4]),[ng^3,ns^3]);
    gNNNx=reshape(permute(reshape(gNNx(:)*gN(:)',[ng^2,ns^2,ng,ns]),[1,3,2,4]),[ng^3,ns^3]);
    gNNNy=reshape(permute(reshape(gNNy(:)*gN(:)',[ng^2,ns^2,ng,ns]),[1,3,2,4]),[ng^3,ns^3]);
    gNNNz=reshape(permute(reshape(gNN(:)*gNx(:)',[ng^2,ns^2,ng,ns]),[1,3,2,4]),[ng^3,ns^3]);
    gfs=permute(cat(3,gNNN,gNNNx,gNNNy,gNNNz),[2,3,1]);
  otherwise
    error('Dimension not implemented.');
end
wgfs=gfs.*repmat(permute(gw,[2,3,1]),[size(gfs,1),dim+1]);
