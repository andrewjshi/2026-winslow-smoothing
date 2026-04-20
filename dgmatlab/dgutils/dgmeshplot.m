function hh=dgmeshplot(msh,opts,expr)
%DGMESHPLOT  Plot DG mesh structure
%
%    Syntax: dgmeshplot(msh,[opts,expr])
%
%    msh:        mesh structure
%    opts:       (logical)
%      opts(1):  plot nodes (default=0)
%      opts(2):  visualize partitions (default=1)
%      opts(3):  plot triangle numbers (default=0)
%    expr:       expression for element inclusion in terms
%                of p (3-D only, default=[])

if nargin<2 | isempty(opts), opts=[0]; end
if length(opts)<2, opts=[opts,isfield(msh,'tpart')]; end
if length(opts)<3, opts=[opts,0]; end
if nargin<3, expr=[]; end

p=msh.p';
t=msh.t'+1;
dim=size(p,2);
if msh.eltype == t_block
    node_plot_order_2d=[1,2,4,3];
else
    node_plot_order_2d=[1,2,3];
end

hh=[];
switch dim
 case 1
  hh=[hh;plot(msh.p,0*msh.p,'.-','linewidth',2,'color','k','markersize',24)];
 case 2
  if opts(2) & ~isempty(msh.tpart)
    ntpart=double(max(msh.tpart(1,:))+1);
    cols=hsv2rgb([[1:ntpart]/ntpart;.5*ones(1,ntpart);ones(1,ntpart)]');
    % Pseudo-random
    jump=ceil(ntpart/2);
    pick=reshape(reshape(1:ceil(ntpart/jump)*jump,jump,[])',[],1);
    pick=pick(pick<=ntpart);
    cols=cols(pick,:);
%    col=cols(msh.tpart(1,:)+1,:);
    col=cols(double(msh.tpart(1,:))+1,:);
    pars={'facecolor','flat','facevertexcdata',col,'edgecolor','k'};
  else
    pars={'facecolor',[.8,1,.8],'edgecolor','k'};
  end
  clf,hh=[hh;patch('faces',t(:,node_plot_order_2d),'vertices',p,pars{:})];
  view(2),axis equal
 case 3
  if ~isempty(expr) & ~ischar(expr)
    t2t=msh.t2t';
    tri1=[];
    %    faces=[2,3,4; 1,4,3; 1,2,4; 1,3,2];
    faces=mkfacemap(3,msh.eltype)';
    for it=1:size(t,1)
      for j=1:size(t2t,2)
        bndnbr=-t2t(it,j);
        if bndnbr>0
          if ismember(bndnbr,expr)
            tri1=[tri1;t(it,faces(j,:))];
          end
        end
      end
    end
  else
    if ~msh.eltype==t_simplex, error('not implemented'); end
    tri1=surftri(p,t);
    if ~isempty(expr)
      incl=find(eval(expr));
      tincl=any(ismember(t,incl),2);
      t1=t(tincl,:);
      tri1=tri1(any(ismember(tri1,incl),2),:);
      tri2=surftri(p,t1);
      tri2=setdiff(tri2,tri1,'rows');
      h=trimesh(tri2(:,node_plot_order_2d),p(:,1),p(:,2),p(:,3));
      hh=[hh;h];
      set(h,'facecolor',[.6,.8,1],'edgecolor','k');
      hold on
    end
  end
  h=trimesh(tri1(:,node_plot_order_2d),p(:,1),p(:,2),p(:,3));
  hh=[hh;h];
  hold off
  set(h,'facecolor',[.8,1,.8],'edgecolor','k');
  axis equal
  
  dragzoom
  cameratoolbar
end

if opts(1)
  if dim==1
    xx=squeeze(msh.p1);
    yy=0*xx;
    zz=0*xx;
  elseif dim==2
    xx=squeeze(msh.p1(:,1,:));
    yy=squeeze(msh.p1(:,2,:));
    zz=0*xx;
  else
    p1=msh.p1;
    if ~isempty(expr)
      if ischar(expr)
        p1=p1(:,:,tincl);
      else
        p2=[];
        t2t=msh.t2t';
        data=dginit(msh);
        for it=1:size(t,1)
          for j=1:size(t2t,2)
            bndnbr=-t2t(it,j);
            if bndnbr>0
              if ismember(bndnbr,expr)
                p2=[p2; p1(data.egix(:,j)+1,:,it)];
              end
            end
          end
        end
        p1=p2;
      end
    end
    xx=squeeze(p1(:,1,:));
    yy=squeeze(p1(:,2,:));
    zz=squeeze(p1(:,3,:));
  end
  line(xx(:),yy(:),zz(:),'lines','n','marker','.','markersi',16,'col','b');
end

if opts(3)
  if dim==2
    for it=1:size(t,1)
      pmid=mean(p(t(it,:),:),1);
      txtpars={'fontname','times','fontsize',12,'horizontala','center'};
      text(pmid(1),pmid(2),num2str(it),txtpars{:});
    end
  elseif dim==3
    t2t=msh.t2t';
    for it=1:size(t,1)
      for j=1:size(t2t,2)
        if ismember(-t2t(it,j),expr)
          faces=mkfacemap(3,msh.eltype);
          pmid=mean(p(t(it,faces(:,j)),:),1);
          txtpars={'fontname','times','fontsize',12,'horizontala','center'};
          text(pmid(1),pmid(2),pmid(3),[num2str(it),',',num2str(j)],txtpars{:});
        end
      end
    end
    set(findobj(gca,'type','patch'),'facealpha',.5)
  end
end

if nargout<1, clear hh; end
