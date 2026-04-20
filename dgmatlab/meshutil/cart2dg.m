function [p,t,p1]=cart2dg(porder,X,Y,Z,opts)

if nargin==3
  opts=struct;
  dim=2;
elseif nargin==4 & isstruct(Z)
  opts=Z;
  dim=2;
elseif nargin==4
  opts=struct;
  dim=3;
elseif nargin==5
  dim=3;
end

if ~isfield(opts,'ind'), opts.ind=0; end
if ~isfield(opts,'eltype'), opts.eltype=t_simplex; end

if dim==2
  [m,n]=size(X);
  XX=permute(cat(3,X,Y),[3,1,2]);
  X1=XX(:,1:porder:m,1:porder:n);
  if opts.eltype==t_simplex
      X1T=block2tets(X1,opts.ind);
  elseif opts.eltype==t_block
      X1T=block2hexes(X1);
  end
  p1=reshape(X1T,2,[])';
  t1=reshape(1:prod(size(X1T))/2,size(X1T,2),size(X1T,3))';
  if opts.eltype==t_simplex
      flip=simpvol(p1,t1)<0;
      t1(flip,[1,2])=t1(flip,[2,1]);
  end
  
  [X0,Y0]=ndgrid(1:m,1:n);
  XX0=permute(cat(3,X0,Y0),[3,1,2]);
  X0=XX0(:,1:porder:m,1:porder:n);
  if opts.eltype==t_simplex
      X0T=block2tets(X0,opts.ind);
  elseif opts.eltype==t_block
      X0T=block2hexes(X0);
  end
  p0=reshape(X0T,2,[])';
elseif dim==3
  [m,n,o]=size(X);

  XX=permute(cat(4,X,Y,Z),[4,1,2,3]);
  X1=XX(:,1:porder:m,1:porder:n,1:porder:o);
  X1T=block2tets(X1,opts.ind);
  p1=reshape(X1T,3,[])';
  t1=reshape(1:prod(size(X1T))/3,4,size(X1T,3))';
  flip=simpvol(p1,t1)<0;
  t1(flip,[1,2])=t1(flip,[2,1]);

  [X0,Y0,Z0]=ndgrid(1:m,1:n,1:o);
  XX0=permute(cat(4,X0,Y0,Z0),[4,1,2,3]);
  X0=XX0(:,1:porder:m,1:porder:n,1:porder:o);
  X0T=block2tets(X0,opts.ind);
  p0=reshape(X0T,3,[])';
end

s=lagrangepnts(porder,dim,opts.eltype);
ns=size(s,1);

nt1=size(t1,1);
dgp0=zeros(ns,dim,nt1);
for comp=1:dim
  for elnode=1:size(s,2)
    dp0=s(:,elnode)*p0(t1(:,elnode),comp)';
    dgp0(:,comp,:)=dgp0(:,comp,:)+permute(dp0,[1,3,2]);
  end
end

dgp0=round(dgp0);
if dim==2
  ix=sub2ind([m,n],dgp0(:,1,:),dgp0(:,2,:));
  dgp(:,1,:)=X(ix);
  dgp(:,2,:)=Y(ix);
elseif dim==3
  ix=sub2ind([m,n,o],dgp0(:,1,:),dgp0(:,2,:),dgp0(:,3,:));
  dgp(:,1,:)=X(ix);
  dgp(:,2,:)=Y(ix);
  dgp(:,3,:)=Z(ix);
end

[p1,t1,pix]=fixmesh(p1,t1,1e-8);

valid_tet = find(simpvol(p1,t1) > 1e-10);
t1 = t1(valid_tet, :);
dgp = dgp(:, :, valid_tet);

p=p1; t=t1; p1=dgp;
