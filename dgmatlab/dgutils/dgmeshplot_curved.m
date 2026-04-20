function hh=dgmeshplot_curved(msh,nref,pltnodes,pltnbrs)

if nargin<2, nref=0; end
if nargin<3, pltnodes=0; end
if nargin<4, pltnbrs=0; end

col=[.8,1,.8];

if ~isfield(msh,'p1'), error('No msh.p1 -- run nodealloc fist.'); end

p1=msh.p1;
t1=msh.tlocal+1;
s1=msh.s;
porder=double(msh.porder);
dim=size(p1,2);
[nf,nt]=size(msh.t);

if dim~=2, error('Only 2-D.'); end

clf,axis equal,axis off
hh=zeros(nt,1);
patchpars={'facecolor',col,'edgecolor','k'};
switch msh.eltype
  case t_simplex
    if nref>0
        A0=pmonomial(s1(:,1:dim),porder);
        [s1,t1]=uniref(s1,t1,nref);
        A=pmonomial(s1(:,1:dim),porder)/A0;
        p1=reshape(A*reshape(p1,size(A0,1),[]),size(A,1),dim,[]);
    end
    e=boundedges(s1(:,1:2),t1);
    e1=segcollect(e);
    for it=1:nt
        px=p1(:,1,it);
        py=p1(:,2,it);
        hh(it)=patch(px(e1{1}'),py(e1{1}'),0.0*e1{1}',patchpars{:});
    end
  case t_block
    ninterp=porder*2^nref+1;
    V=plegendre(linspace(-1,1,ninterp)',porder)/plegendre(2*msh.s0-1,porder);
    for it=1:nt
        x=reshape(msh.p1(:,1,it),[porder+1,porder+1]);
        y=reshape(msh.p1(:,2,it),[porder+1,porder+1]);
        xc=V*[x(:,1),x(end,:)',x(end:-1:1,end),x(1,end:-1:1)'];
        yc=V*[y(:,1),y(end,:)',y(end:-1:1,end),y(1,end:-1:1)'];
        hh(it)=patch(xc(:),yc(:),0*xc(:),patchpars{:});
    end
  otherwise
    error('Unknown element type');
end

if pltnodes
  x1=msh.p1(:,1,:);
  y1=msh.p1(:,2,:);
  hh1=line(x1(:),y1(:),'linestyle','none','marker','.','col','b');
  hh=[hh; hh1];
end

if pltnbrs==1
  for it=1:size(msh.t,2)
    pmid=mean(msh.p(:,msh.t(:,it)+1),2);
    hh1=text(pmid(1),pmid(2),int2str(it),'horiz','center','verti','middle');
    hh=[hh; hh1];
  end
end

if nargout<1, clear hh; end
