function p=quadsmooth(p,q,bndix,n)

if nargin<4, n=1; end

edges=[q(:,1),q(:,2); q(:,2),q(:,3); q(:,3),q(:,4); q(:,4),q(:,1)];
%edges=[q(:,1),q(:,2); q(:,2),q(:,3); q(:,3),q(:,4); q(:,4),q(:,1); q(:,1),q(:,3); q(:,2),q(:,4)];
edges=unique(sort(edges,2),'rows');

for i=1:n
  px=p(:,1);
  py=p(:,2);
  
  nnew=full(sparse(edges,1,1));
  pxnew=full(sparse(edges,1,px(fliplr(edges))))./nnew;
  pynew=full(sparse(edges,1,py(fliplr(edges))))./nnew;
  
  pxnew(bndix)=px(bndix);
  pynew(bndix)=py(bndix);
  
  p=[pxnew,pynew];
  % clf,meshplotcol(p,q,quadqual(p,q)),pause
end
