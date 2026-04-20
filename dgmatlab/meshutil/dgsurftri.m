function [p,tri] = dgsurftri(msh, nref, bnds, fixmsh)
%DGSURFTRI  Return surface mesh for given boundaries
%   [P,TRI] = DGSURFTRI(MSH, NREF, BNDS)
%
%      NREF = -1: Straight-sided mesh
%           =  0: Interpolate DG nodes
%           >  0: Refine NREF times (on true curved geometry)
%
%      BNDS = list of boundary numbers (default all)
    
if nargin<3 | isempty(bnds), bnds=1:(-min(msh.t2t(:))); end
if nargin<4, fixmsh=true; end

nt=size(msh.t,2);
porder=double(msh.porder);
eltype=t_simplex;

if nref == -1  % Straight sided mesh
    p=msh.p';
    t=double(msh.t')+1;
    t2t=double(msh.t2t');
    tri1=[];
    faces=mkfacemap(3,msh.eltype)';
    for it=1:size(t,1)
        for j=1:size(t2t,2)
            bndnbr=-t2t(it,j);
            if bndnbr>0
                if ismember(bndnbr,bnds)
                    tri1=[tri1;t(it,faces(j,:))];
                end
            end
        end
    end
    [p,tri] = fixmesh(p, tri1);
    return;
end

if porder==0
  [ss,tt]=lagrangepnts(1,2);
  pp=zeros(dim+1,dim,nt);
  for i=1:dim
    pp(:,i,:)=reshape(msh.p(i,msh.t+1),4,1,nt);
  end
else
  ss=msh.sbnd;
  tt=msh.tbndlocal;
  if min(tt(:)) == 0
      tt = tt + 1;
  end
  pp=msh.p1;
end
ns=size(ss,1);

if nargin<2 | isempty(nref), nref = 0; end

if msh.eltype == t_block
    node_plot_order_2d=[1,2,4,3];
else
    node_plot_order_2d=[1,2,3];
end

data=dginit(msh);

switch msh.eltype
  case t_simplex
    tbnd=msh.tbndlocal;
    if min(tbnd(:)) == 0
        tbnd = tbnd + 1;
    end
    if nref>0
        A0=pmonomial(ss(:,1:2),porder);
        [ss,tbnd]=lagrangepnts(porder*2^nref,2);
        A=pmonomial(ss(:,1:2),porder)/A0;
    end
    tbnd=double(tbnd);
    nbnd=size(tbnd,1);
    ns=size(ss,1);
    
    nplt=sum(ismember(-msh.t2t(:),bnds));
    if nplt==0, hh=[]; return; end
    
    nsbnd=size(msh.sbnd,1);
    uplt=zeros(nsbnd,nplt);
    pplt=zeros(nsbnd,3,nplt);
    cix=0;
    bndsvec = accumarray(bnds(:), 1, [-min(msh.t2t(:)),1]);
    for it=1:nt
        for j=1:4
            bndnbr=-msh.t2t(j,it);
            if bndnbr>=1 & bndsvec(bndnbr)
                cegix=data.egix(:,j)+1;
                cp=msh.p1(cegix,:,it);
                cix=cix+1;
                pplt(:,:,cix)=cp;
            end
        end
    end
    pplt=reshape(pplt,nsbnd,3*nplt);
    if nref>0
        pplt=A*pplt;
    end
    ntri=size(pplt,2) / 3;
    tplt=kron(ones(ntri,1),tbnd)+kron((0:ntri-1)'*ns,ones(nbnd,3));

    pplt=reshape(permute(reshape(pplt,ns,3,ntri),[1,3,2]),[],3);
  case t_block
    tbnd=msh.tbndlocal;
    if min(tbnd(:)) == 0
        tbnd = tbnd + 1;
    end
    if nref>0
        error('not implemented');
    end
    tbnd=double(tbnd);
    nbnd=size(tbnd,1);
    ns=size(ss,1);
    
    nplt=sum(ismember(-msh.t2t(:),bnds));
    if nplt==0, hh=[]; return; end
    
    nsbnd=size(msh.sbnd,1);
    pplt=zeros(nsbnd,3,nplt);
    cix=0;
    for it=1:nt
        for j=1:6
            bndnbr=-msh.t2t(j,it);
            if ismember(bndnbr,bnds)
                cegix=data.egix(:,j)+1;
                cp=msh.p1(cegix,:,it);
                cix=cix+1;
                pplt(:,:,cix)=cp;
            end
        end
    end
    pplt=reshape(pplt,nsbnd,3*nplt);
    if nref>0
        error('not implemented');
    end
    nquad=size(pplt,2);
    tplt=kron(ones(nquad,1),tbnd)+kron((0:nquad-1)'*ns,ones(nbnd,4));
    pplt=reshape(permute(reshape(pplt,ns,3,nquad),[1,3,2]),[],3);
    tplt=tplt(:,node_plot_order_2d);
end
 
if fixmsh
    [p,tri] = fixmesh(pplt, tplt);
else
    p = pplt; tri = tplt;
end
