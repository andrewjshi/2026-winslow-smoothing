function msh=nodealloc(msh,porder)

if isscalar(porder)
    porder=double(porder);
    if msh.eltype == t_simplex
        msh.s0 = (0:porder)'/porder;
    elseif msh.eltype == t_block
        msh.s0=(lobattopoints(porder+1)+1)/2;
    else
        error('Assertion: Unknown element type');
    end
else % Assume porder is s0
    msh.s0 = porder;
    porder = numel(msh.s0) - 1;
end
msh.porder=int32(porder);

dim=size(msh.p,1);
[msh.s,msh.tlocal,msh.sbnd,msh.tbndlocal]=lagrangepnts(porder,dim,msh.eltype,msh.s0);
msh.tlocal = msh.tlocal - 1;
msh.tbndlocal = msh.tbndlocal - 1;
[nf,nt]=size(msh.t2t);
ns=size(msh.s,1);

if ~isfield(msh,'ecurved')
  error('No ecurved, run mshcurved.');
end
if ~isfield(msh,'ncurved')
  error('No ncurved, run mshcurved.');
end
if ~isfield(msh,'tcurved')
  msh.tcurved=any(msh.ecurved,1);
end

% Constraint on having curved elements first is removed for now
% if any(diff(msh.tcurved)>0)
%   error('All curved elements must be located first in msh.t, run mshfix.');
% end

% Allocate nodes
msh.p1=zeros(ns,dim,nt);
for comp=1:dim
  for elnode=1:size(msh.s,2)
    dp=msh.s(:,elnode)*msh.p(comp,double(msh.t(elnode,:))+1);
    msh.p1(:,comp,:)=msh.p1(:,comp,:)+permute(dp,[1,3,2]);
  end
end

% If got distance function, project and smooth curved boundary elements
if isfield(msh,'fd') & ~isempty(msh.fd)
  msh=smoothdgnodes(msh);
end
