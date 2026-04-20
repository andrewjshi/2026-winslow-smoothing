function hh=dgplot1d(msh,u,utype,clim,nref,pltmesh)

nt=size(msh.t,2);
ns=size(msh.s,1);

porder=double(msh.porder);

if porder==0
  [ss,tt]=lagrangepnts(1,1);
  pp=zeros(2,1,nt);
  pp(:,1,:)=reshape(msh.p(1,msh.t+1),2,1,nt);
  u=repmat(u,[2,1,1]);
else
  ss=msh.s;
  tt=msh.tlocal;
  if min(tt(:)) == 0
      tt = tt + 1;
  end
  pp=msh.p1;
end

if size(u,1)==1
  u=repmat(u,[ns,1,1]);
end

if nargin<5 | isempty(nref), nref=ceil(log2(max(porder,1))); end

if nref>0
  A0=pmonomial(ss(:,1),porder);
  [ss,tt]=uniref(ss,tt,nref);
  A=pmonomial(ss(:,1),porder)/A0;

  nss=size(ss,1);
  sz=size(pp); if length(sz)==2, sz=[sz,1]; end
  pp=reshape(A*reshape(pp,ns,sz(2)*sz(3)),[nss,sz(2),sz(3)]);
  sz=size(u); if length(sz)==2, sz=[sz,1]; end
  u=reshape(A*reshape(u,ns,sz(2)*sz(3)),[nss,sz(2),sz(3)]);
end

nss=size(ss,1);
p1vis=reshape(permute(pp,[1,3,2]),[nss*nt,1]);
tvis=kron(ones(nt,1,'int32'),tt)+kron(nss*int32(0:nt-1)',0*tt+1);

if nargin>=3 & ~isempty(utype)
  u=dgeval(u,utype,1);
else
  u=u(:,1,:);
end

if strcmp(get(gca,'nextplot'),'replace')
    cla
end

if nargin>=4 & ~isempty(clim)
else
  minu=min(u(:));
  maxu=max(u(:));
  du=maxu-minu;
  if du==0
    minu=0; maxu=1; du=1;
  end
  minu=minu-.1*du;
  maxu=maxu+.1*du;
  clim=[minu,maxu];
end
set(gca,'ylim',clim);

if nargin>=6 & ~isempty(pltmesh) & pltmesh
  line([msh.p;msh.p],[0*msh.p+clim(1); 0*msh.p+clim(2)], ...
       'color',[.9,.9,.9])
end

hh=line(p1vis(tvis'),u(tvis'),'color','k');

set(gcf,'rend','z');
drawnow
