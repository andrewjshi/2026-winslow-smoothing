function t2perm=mkt2perm(t,t2t,t2n,eltype)

if nargin<4, eltype=t_simplex; end

[nt,nf]=size(t2t);
dim=nv2dim(size(t,2),eltype);
t2perm=ones(size(t2t),'int32');

if dim==1 || dim==2
  return;
elseif dim~=3
  error('Dimension not implemented');
end

map=mkfacemap(dim,eltype);
for it1=1:nt
  for j=1:nf
    it2=t2t(it1,j);
    if it2>=1
      k=mod(t2n(it1,j),16);
      f1=t(it1,map(:,j));
      f2=t(it2,map(:,k));
      t2perm(it1,j)=int32(find(f2==f1(1)));
    end
  end
end
