function [t2t,t2n]=mkt2t(t,perminfo,eltype)

if nargin<3, eltype=t_simplex; end

nt=size(t,1);
edges=mkfaces(t,eltype);
dim=nv2dim(size(t,2),eltype);
nf=size(mkfacemap(dim,eltype),2);
ts=[repmat(int32(1:nt),1,nf); kron(int32(1:nf),ones(1,nt,'int32'))]';

edges=sort(edges,2);
[foo,foo,jx]=unique(edges,'rows');

[jx,ix]=sort(jx);
ts=ts(ix,:);

ix=find(diff(jx)==0);
ts1=ts(ix,:);
ts2=ts(ix+1,:);

t2t=zeros(nt,nf,'int32');
t2t(ts1(:,1)+nt*(ts1(:,2)-1))=ts2(:,1);
t2t(ts2(:,1)+nt*(ts2(:,2)-1))=ts1(:,1);

if nargout>=2
  t2n=zeros(nt,nf,'int32');
  t2n(ts1(:,1)+nt*(ts1(:,2)-1))=ts2(:,2);
  t2n(ts2(:,1)+nt*(ts2(:,2)-1))=ts1(:,2);  
  
  if nargin>=2 & perminfo
    t2perm=mkt2perm(t,t2t,t2n,eltype);
    t2n=t2n+2^4*(t2perm-1)+2^7;
  end
end
