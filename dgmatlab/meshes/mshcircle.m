function msh=mshcircle(nref)

if nargin<1, nref=0; end

fd=inline('sqrt(sum(p.^2,2))-1','p');
fdargs={};

phi=2*pi*(0:7)'/8;
p=[0,0;cos(phi),sin(phi)];
t=[1,2,3;1,3,4;1,4,5;1,5,6;1,6,7;1,7,8;1,8,9;1,9,2];

for iref=1:nref
  [p,t]=uniref(p,t,1);
  p=bndproj(p,t,fd,fdargs{:});
end

bndexpr={'true | p(:,1)'};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,1);
msh=mshfix(msh);
