function [msh,perm]=mshreorder(msh, order, varargin)

if isnumeric(order)
  perm=order;
elseif isa(order, 'function_handle')
  A=dualadjacency(msh.t2t'+1);
  perm=feval(order,A);
elseif strcmp(order, 'weight')
    if nargin >= 3
        bndcnds = varargin{1};
    else
        bndcnds = [2,1];
    end
    phys=physinit(msh, bndcnds, [1.4,.1], [1e3,.72]);
    dt=1e-3;
    data=dginit(msh);
    u=freestream(msh,phys);
    perm=matweightperm(u, msh, data, phys, @dgnavierstokes, dt, [], 0);
else
    error('Invalid order');
end

msh.t=msh.t(:,perm);
msh.t2n=msh.t2n(:,perm);

msh.t2t=msh.t2t(:,perm);
ix=msh.t2t>=0;
[foo,iperm]=sort(perm);
msh.t2t(ix)=iperm(msh.t2t(ix)+1)-1;

if isfield(msh,'ecurved')
  msh.ecurved=msh.ecurved(:,perm);
end
if isfield(msh,'tcurved')
  msh.tcurved=msh.tcurved(:,perm);
end

if isfield(msh,'p1')
  msh.p1=msh.p1(:,:,perm);
end
if isfield(msh,'dd')
  msh.dd=msh.dd(:,:,perm);
end
