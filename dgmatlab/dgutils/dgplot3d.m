function hh=dgplot3d(msh,u,utype,clim,nref,pltmesh,bnds)

if nargin<7 | isempty(bnds), bnds=1:(-min(msh.t2t(:))); end
if nargin<6 | isempty(pltmesh), pltmesh=1; end

nt=size(msh.t,2);
porder=double(msh.porder);
eltype=t_simplex;

if porder==0
  [ss,tt]=lagrangepnts(1,2);
  pp=zeros(dim+1,dim,nt);
  for i=1:dim
    pp(:,i,:)=reshape(msh.p(i,msh.t+1),4,1,nt);
  end
  u=repmat(u,[4,1,1]);
else
  ss=msh.sbnd;
  tt=msh.tbndlocal;
  if min(tt(:)) == 0
      tt = tt + 1;
  end
  pp=msh.p1;
end
ns=size(ss,1);

if size(u,1)==1
  u=repmat(u,[ns,1,1]);
end

if nargin<5 | isempty(nref), nref=ceil(log2(max(porder,1))); end

if nargin>=3 & ~isempty(utype)
  u=dgeval(u,utype,3);
else
  u=u(:,1,:);
end

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
                cu=u(cegix,1,it);
                cp=msh.p1(cegix,:,it);
                cix=cix+1;
                uplt(:,cix)=cu;
                pplt(:,:,cix)=cp;
                %      uplt=[uplt, cu];
                %      pplt=[pplt, cp];
            end
        end
    end
    pplt=reshape(pplt,nsbnd,3*nplt);
    if nref>0
        uplt=A*uplt;
        pplt=A*pplt;
    end
    ntri=size(uplt,2);
    tplt=kron(ones(ntri,1),tbnd)+kron((0:ntri-1)'*ns,ones(nbnd,3));

    pplt=reshape(permute(reshape(pplt,ns,3,ntri),[1,3,2]),[],3);
    uplt=uplt(:);

    hh=[];
    if ~ishold, cla; end
    h=patch('vertices',pplt,'faces',tplt,'cdata',uplt, ...
            'facecol','interp','edgec','none');
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
    uplt=zeros(nsbnd,nplt);
    pplt=zeros(nsbnd,3,nplt);
    cix=0;
    for it=1:nt
        for j=1:6
            bndnbr=-msh.t2t(j,it);
            if ismember(bndnbr,bnds)
                cegix=data.egix(:,j)+1;
                cu=u(cegix,1,it);
                cp=msh.p1(cegix,:,it);
                cix=cix+1;
                uplt(:,cix)=cu;
                pplt(:,:,cix)=cp;
                %      uplt=[uplt, cu];
                %      pplt=[pplt, cp];
            end
        end
    end
    pplt=reshape(pplt,nsbnd,3*nplt);
    if nref>0
        error('not implemented');
    end
    nquad=size(uplt,2);
    tplt=kron(ones(nquad,1),tbnd)+kron((0:nquad-1)'*ns,ones(nbnd,4));
    pplt=reshape(permute(reshape(pplt,ns,3,nquad),[1,3,2]),[],3);
    uplt=uplt(:);
    tplt=tplt(:,node_plot_order_2d);
    
    hh=[];
    if ~ishold, cla; end
    h=patch('vertices',pplt,'faces',tplt,'cdata',uplt, ...
            'facecol','interp','edgec','none');
end
    
hh=[hh;h];
if nargin>=4 & ~isempty(clim)
  if isequal(clim,1)
    clim=[min(u(:)),max(u(:))];
  else
    set(gca,'clim',clim);
  end
end

if pltmesh==1 % Straight sided mesh
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
  h=patch('vertices',p,'faces',tri1(:,node_plot_order_2d),'facecol','none','edgec',[0,0,0]);
  hh=[hh;h];
elseif pltmesh==2 % Curved mesh
  e=boundedges(pplt,tplt,eltype);
  ppltx=pplt(:,1);
  pplty=pplt(:,2);
  ppltz=pplt(:,3);
%  h=line(ppltx(e'),pplty(e'),ppltz(e'),'color',[0,0,0]);
  h=patch('vertices',pplt,'faces',e,'facecol','none','edgec',[0,0,0]);
  hh=[hh;h];
end

axis equal
%view(3)

dragzoom
cameratoolbar
  
set(gcf,'rend','z');
colorbar,drawnow
