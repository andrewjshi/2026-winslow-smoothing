function msh=mshrae1(porder,resolution)

load rae_mesh
X=rae(1:resolution:end,1:resolution:end,:);
X=permute(X,[3,1,2]);
XT=block2tets(X);

p0=reshape(XT,2,[])';
t0=reshape(1:prod(size(XT))/2,2+1,size(XT,3))';
[p,t]=fixmesh(p0,t0);

bndexpr={'all(sqrt(sum(p.^2,2))<2)','true'};
msh=ml2msh(p,t,bndexpr);

msh.ncurved=logical(ones(1,size(msh.p,2)));
msh.ecurved=logical(ones(3,size(msh.t,2)));
msh.tcurved=logical(ones(1,size(msh.t,2)));

[msh.s,msh.tlocal]=lagrangepnts(porder,2);
[msh.sbnd,msh.tbndlocal]=lagrangepnts(porder,1);

error('No DG nodes created.');

% To Do: Create msh.p1 using nodes in rae
% ...
ns=size(msh.s,1); nt=size(msh.t,2);
p1=zeros(ns,2,nt);
for it=1:nt
  
end
