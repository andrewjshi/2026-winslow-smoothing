function h2=tetplotcurved(msh,bnd,ix,color, nref)

if nargin<3, ix=find(any(ismember(msh.t2t,-bnd),1)); end
if nargin<4, color = [1, .5, .5]; end
if nargin<5 | isempty(nref), nref=2; end

%h1=dgplot(msh,msh.p1*0,1,[],nref,2,bnd);
%set(findobj(h1,'type','patch'),'facecol',[.8,1,.8])

if ~isempty(ix)
    p = msh.p';
    t = msh.t'+1;
    t = t(ix,:);
    t2t = mkt2t(t);
    tri1 = surftri(p,t);
    
    msh1=msh;
    msh1.t=msh1.t(:,ix);
    msh1.p1=msh1.p1(:,:,ix);
    msh1.t2t = int32(t2t'-1);
  
    hold on
    view(-50, 10)
    h2=dgplot(msh1,msh1.p1*0,1,[],nref,2,1);
    hold off
    %set(findobj(h2,'type','patch'),'facecol',[1,.5,.5])
    set(findobj(h2,'type','patch'),'facecol',color)
    %set(findobj(h2,'type','patch'),'facealpha',.1)
end
