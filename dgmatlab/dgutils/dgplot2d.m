function dgplot2d(msh,u,utype,clim,nref,pltmesh)

if nargin<3, utype=1; end
if nargin<4, clim=[]; end
if nargin<5, nref=0; end
if nargin<6, pltmesh=0; end

[nv,nt]=size(msh.t);
dim=size(msh.p,1);
ns=size(msh.s,1);
porder=double(msh.porder);

if porder==0
  [ss,tt]=lagrangepnts(1,dim,msh.eltype);
  pp=zeros(nv,2,nt);
  pp(:,1,:)=reshape(msh.p(1,msh.t+1),nv,1,nt);
  pp(:,2,:)=reshape(msh.p(2,msh.t+1),nv,1,nt);
  u=repmat(u,[nv,1,1]);
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
u=dgeval(u,utype,dim);

cla
patchpars={'facecol','interp','edgec','none'};
meshcol=[0,0,0];
switch msh.eltype
  case t_simplex
    if nref>0
        A0=pmonomial(ss(:,1:2),porder);
        [ss,tt]=uniref(ss,tt,nref);
        A=pmonomial(ss(:,1:2),porder)/A0;
        
        nss=size(ss,1);
        sz=size(pp); if length(sz)==2, sz=[sz,1]; end
        pp=reshape(A*reshape(pp,ns,sz(2)*sz(3)),[nss,sz(2),sz(3)]);
        sz=size(u); if length(sz)==2, sz=[sz,1]; end
        u=reshape(A*reshape(u,ns,sz(2)*sz(3)),[nss,sz(2),sz(3)]);
    end
    
    nss=size(ss,1);
    p1vis=reshape(permute(pp,[1,3,2]),[nss*nt,2]);
    tvis=kron(ones(nt,1,'int32'),tt)+kron(nss*int32(0:nt-1)',0*tt+1);
    
    patch('vertices',p1vis,'faces',tvis,'cdata',squeeze(u),patchpars{:});
    
    switch pltmesh
      case 0 % No mesh
      case 1 % Plot straight-sided mesh
        patch('vertices',msh.p','faces',double(msh.t)'+1,'facecolor','none','edgecolor',meshcol);
      case 2 % Plot curved mesh
        e=boundedges(p1vis,tvis);
        ppltx=p1vis(:,1);
        pplty=p1vis(:,2);
        line(ppltx(e'),pplty(e'),'color',meshcol);
      case 3 % Plot the refined mesh
        patch('vertices',p1vis,'faces',tvis,'facecol','none','edgecolor',meshcol);
      otherwise
        error('Unknown mesh plot type.');
    end
  case t_block
    ninterp=porder*2^nref+1;
    V=plegendre(linspace(-1,1,ninterp)',porder)/plegendre(2*msh.s0-1,porder);
    for it=1:size(msh.t,2)
        x=V*reshape(msh.p1(:,1,it),[porder+1,porder+1])*V';
        y=V*reshape(msh.p1(:,2,it),[porder+1,porder+1])*V';
        c=V*reshape(u(:,:,it),[porder+1,porder+1])*V';
        surface(x,y,0*x,c,patchpars{:});
        switch pltmesh
          case 0 % No mesh
          case 1 % Plot straight-sided mesh
            patch('vertices',msh.p','faces',double(msh.t([1,2,4,3],:))'+1,'facecolor','none','edgecolor',meshcol);
          case 2 % Plot curved mesh
            xc=[x(:,1);x(end,:)';x(end:-1:1,end);x(1,end:-1:1)'];
            yc=[y(:,1);y(end,:)';y(end:-1:1,end);y(1,end:-1:1)'];
            line(xc,yc,'color',[0,0,0]);
          otherwise
            error('Unknown mesh plot type.');
        end
        if pltmesh==2
            % Plot curved mesh
        end
    end
  otherwise
    error('Unknown element type');
end
    
if ~isempty(clim)
    set(gca,'clim',clim);
end
set(gcf,'rend','z');
%colorbar,drawnow
%drawnow
