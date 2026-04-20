function msh=mshchannel(m,n,r,t)

if numel(t)==1, t=[t,t]; end
if numel(r)==1, r=[r,r]; end

XT=block2tets(block([m,n],[t,r],[0,0;1,0;0,.2;1,.2]));
p=reshape(XT,2,[])';
t=reshape(1:numel(XT)/2,3,size(XT,3))';
[p,t]=fixmesh(p,t);

bndexpr={'all(p(:,2)<1e-6)','all(p(:,1)>1-1e-6)', ...
         'all(p(:,2)>.2-1e-6)','all(p(:,1)<1e-6)'};
msh=ml2msh(p,t,bndexpr);
msh=mshcurved(msh,[]);
