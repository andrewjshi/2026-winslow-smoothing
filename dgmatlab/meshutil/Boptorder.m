function pp=Boptorder(X0)

[X,ix,jx]=unique(X0,'rows');
%p=minflip_greedy(X);
p=minflip_twoopt(X);
pp=[];
for ii=1:length(p)
  pp=[pp;find(jx==p(ii))];
end

function p=minflip_greedy(X)

np=size(X,2);
c=sum(X,2);

cc=zeros(1,np);
for i=1:np
  cc(i)=sum(~~sum(X(X(:,i)~=0,:),1));
end

[foo,ix]=min(cc);
ixx=find(X(:,ix)~=0);
[foo,ixxx]=min(c(ixx));
p=ixx(ixxx);

XX=[X,(1:size(X,1))'];
cX=X(p,1:np);
XX(p,:)=[];

while size(XX,1)>0
  [foo,minix]=min(sum(repmat(cX,size(XX,1),1)~=XX(:,1:np),2));
  p(end+1)=XX(minix,np+1);
  cX=XX(minix,1:np);
  XX(minix,:)=[];
end

function p=minflip_twoopt(X)

n=size(X,1);
W=zeros(n+1,n+1);

for i=1:n
  for j=1:n
    W(i,j)=sum(diff(X([i,j],:),[],1)==1);
  end
  W(n+1,i)=sum(X(i,:));
end

if n>0
  p0=minflip_greedy(X);
  p0=[p0,n+1];
  p=twoopt(int32(W),int32(p0));
else
  % Need to do some debugging of this
  W=W+1e6*eye(size(W));
  p=babtsp(int32(W));
end

ix=find(p==n+1);
p=p([ix+1:end,1:ix-1]);
