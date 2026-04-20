function msh=mshsphere(n)

if nargin<1, n=1; end

fd=@dsphere;
fdargs={0,0,0,1};

% To do: Fix

[p,e,t]=tetrahedronmesh(n+1);
p=p.*repmat(sum(p,2)./(sqrt(sum(p.^2,2))+(all(p==0,2))),1,3);

x=p(:,1); y=p(:,2); z=p(:,3); np=size(p,1);
p=[ x, y, z;
   -x, y, z;
    x,-y, z;
   -x,-y, z;
    x, y,-z;
   -x, y,-z;
    x,-y,-z;
   -x,-y,-z];
t=[t+0*np;
   t+1*np;
   t+2*np;
   t+3*np;
   t+4*np;
   t+5*np;
   t+6*np;
   t+7*np];
[p,t]=fixmesh(p,t);

bndexpr={'true | p(:,1)'};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,1);
msh=mshfix(msh);
