function msh=ng2msh(fname,porder)
%NG2MSH NETGEN mesh to 3DG msh
%   Currently assuming 3-D second-order elements

map=[6,3,10,1,5,9,4,8,2,7];
imap(map)=1:10;

[p,t,e,tnbr,enbr]=ng2pt(fname);
t1=t(:,[2,1,3,4]);

[pix,ix1,jx1]=unique(t1);
t1=reshape(jx1,size(t1));
p1=p(pix,:);

msh=ml2msh(p1,t1);
msh=mshcurved(msh,'all');
msh=nodealloc(msh,2);

p1ix=t(:,imap)';
msh.p1=permute(reshape(p(p1ix(:),:),[10,size(t,1),3]),[1,3,2]);

faces=reshape(permute(cat(3,t1(:,[2,3,4]),t1(:,[3,4,1]),t1(:,[4,1,2]),t1(:,[1,2,3])),[3,1,2]),[],3);
[foo,ix,iy]=intersect(sort(faces,2),sort(e(:,1:3),2),'rows');
msh.t2t(ix)=-int32(enbr(iy));

if nargin>=2 & ~isempty(porder) & porder~=2
  p1_2=msh.p1;
  msh=nodealloc(msh,porder);
  msh.p1=dginterp(p1_2,2,porder,3);
end
