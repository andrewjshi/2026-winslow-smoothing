function msh=cart2msh(porder,X,Y,Z,bndexpr,opts)
%CART2MSH
%   2-D:
%     MSH=CART2MSH(PORDER,X,Y,BNDEXPR)
%   3-D:
%     MSH=CART2MSH(PORDER,X,Y,Z,BNDEXPR)
%   General (see CART2DG for creating P,T,P1):
%     MSH=CART2MSH(PORDER,P,T,P1,BNDEXPR)

if nargin<4 | iscell(Z)
  dim=2;
  if nargin>=5
      opts=bndexpr;
  end
  if nargin>=4
      bndexpr=Z;
  end
  ptform=false;
elseif isequal(size(X),size(Y))
  dim=3;
  ptform=false;
else
  p=X; t=Y; p1=Z;
  dim=size(p,2);
  ptform=true;
end

if ~exist('bndexpr'), bndexpr={}; end
if ~exist('opts'), opts=struct; end

if ~ptform
  if dim==2
    [p,t,p1]=cart2dg(porder,X,Y,opts);
  elseif dim==3
    [p,t,p1]=cart2dg(porder,X,Y,Z,opts);
  end
end

msh=ml2msh(p,t,bndexpr,[],[],[],'all');
msh=nodealloc(msh,porder);
msh.p1=p1;

J=geojac(msh);
if any(J(:)<0), warning('  Inverted elements'); end
