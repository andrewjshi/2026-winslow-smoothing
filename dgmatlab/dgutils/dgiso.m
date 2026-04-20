function hh=dgiso(msh,u,utype,up,uptype,lev,nref,cont,normals,chunk)

if nargin<8, cont=0; end
if nargin<9, normals=0; end
if nargin<10, chunk=inf; end

if cont
  u=dg2continuous(msh,u);
  up=dg2continuous(msh,up);
end

if ~isinf(chunk)
  ct=1;
  nt=size(msh.t,2);
  hh=[];
  while ct<=nt
    lt=min(ct+chunk-1,nt);
    ix=ct:lt;
    msh1=msh;
    msh1.t=msh1.t(:,ix);
    msh1.p1=msh1.p1(:,:,ix);
    u1=u(:,:,ix);
    up1=up(:,:,ix);
    hh=[hh,dgiso(msh1,u1,utype,up1,uptype,lev,nref,false,normals,inf)];
    ct=lt+1;
  end
  return;
end

[nv,nt]=size(msh.t);
nt=size(msh.t,2);
ns=size(msh.s,1);
dim=3;

porder=double(msh.porder);

if porder==0
  [ss,tt]=lagrangepnts(1,dim,msh.eltype);
  pp=zeros(nv,dim,nt);
  for i=1:dim
    pp(:,i,:)=reshape(msh.p(i,msh.t+1),nv,1,nt);
  end
  u=repmat(u,[nv,1,1]);
  up=repmat(up,[nv,1,1]);
else
  ss=msh.s;
  tt=msh.tlocal+1;
  if min(tt(:)) == 0
      tt = tt + 1;
  end
  pp=msh.p1;
end

if size(u,1)==1
  u=repmat(u,[ns,1,1]);
  up=repmat(up,[ns,1,1]);
end

if nargin>=3 & ~isempty(utype)
  u=dgeval(u,utype,3);
else
  u=u(:,1,:);
end

if nargin>=5 & ~isempty(uptype)
  up=dgeval(up,uptype,3);
else
  up=up(:,1,:);
end

if nargin<7 | isempty(nref), nref=0; end

switch msh.eltype
  case t_simplex
    if nref>0
        A0=pmonomial(ss(:,1:3),porder);
        [ss,tt]=lagrangepnts(porder*2^nref,3);
        A=pmonomial(ss(:,1:3),porder)/A0;
        
        nss=size(ss,1);
        sz=size(pp); if length(sz)==2, sz=[sz,1]; end
        pp=reshape(A*reshape(pp,ns,sz(2)*sz(3)),[nss,sz(2),sz(3)]);
        sz=size(u); if length(sz)==2, sz=[sz,1]; end
        u=reshape(A*reshape(u,ns,sz(2)*sz(3)),[nss,sz(2),sz(3)]);
        sz=size(up); if length(sz)==2, sz=[sz,1]; end
        up=reshape(A*reshape(up,ns,sz(2)*sz(3)),[nss,sz(2),sz(3)]);
    end

  case t_block
    ninterp=porder*2^nref+1;
    news = linspace(0, 1, ninterp);
    pp = dginterp(msh.p1, msh.s0, news, 3, t_block);
    u = dginterp(u, msh.s0, news, 3, t_block);
    up = dginterp(up, msh.s0, news, 3, t_block);

    msh0 = qmshcube(1,1,1);
    msh0 = nodealloc(msh0, ninterp-1);
    ss = msh0.s;
    tt = msh0.tlocal;
    
    if min(tt(:)) == 0
        tt = tt + 1;
    end
  otherwise
    error('Unknown element type');
end

nss=size(ss,1);
p1vis=reshape(permute(pp,[1,3,2]),[nss*nt,3]);
tvis=kron(ones(nt,1,'int32'),tt)+kron(nss*int32(0:nt-1)',0*tt+1);
    
if msh.eltype == t_block
    tvis = [tvis(:,[5,6,3,7]);
            tvis(:,[5,6,2,3]);
            tvis(:,[5,1,3,2]);
            tvis(:,[6,8,4,7]);
            tvis(:,[6,7,4,3]);
            tvis(:,[6,2,3,4])];
end

u=u(:);
up=up(:);

[X,Y,Z,D]=ptiso(p1vis',int32(tvis'-1),u(tvis),lev,up(tvis));
if isempty(X), hh=[]; return; end

p=[X(:),Y(:),Z(:)];
D=D(:);
t=reshape(1:prod(size(X)),3,[])';

[p,t,pix]=fixmesh(p,t,1e-5);
D=D(pix);

if normals
  nn=ptnormals(p,t);
  hh=patch('vertices',p,'faces',t,'cdata',D,'facecol','interp','edgecol','none', ...
           'vertexnormals',nn);
else
  hh=patch('vertices',p,'faces',t,'cdata',D,'facecol','interp','edgecol','none');
end

function vn=ptnormals(p,t)

t12=p(t(:,2),:)-p(t(:,1),:);
t13=p(t(:,3),:)-p(t(:,1),:);
n=cross(t12,t13);
n=n./repmat(sqrt(sum(n.^2,2)),1,3);

np=size(p,1);
vnx=full(sparse(t(:),1,repmat(n(:,1),3,1),np,1));
vny=full(sparse(t(:),1,repmat(n(:,2),3,1),np,1));
vnz=full(sparse(t(:),1,repmat(n(:,3),3,1),np,1));
vn=[vnx,vny,vnz];
vnnorm=sqrt(sum(vn.^2,2));
vnnorm(vnnorm==0)=1;
vn=vn./repmat(vnnorm,1,3);
vn(:,[2,3])=vn(:,[3,2]);
