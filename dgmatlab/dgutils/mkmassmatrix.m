function M=mkmassmatrix(msh,data,N,block,doinv)

if nargin<4, block=false; end
if nargin<5, doinv=false; end

MD=dgmass(msh,data);
if doinv
  sz=repmat(sqrt(size(MD,1)),1,2);
  for i=1:size(MD,2)
    MD(:,i)=reshape(inv(reshape(MD(:,i),sz)),[],1);
  end
end

if block
  ns=sqrt(size(MD,1));
  nt=size(MD,2);
  MD=reshape(MD,ns,ns,nt);
  M=zeros(N*ns,N*ns,nt);
  for i=1:N
    ix=(1:ns)+ns*(i-1);
    M(ix,ix,:)=MD;
  end
else
  [Dii,Djj]=dgindices(msh,data,'mass',N);
  M=sparse(Dii,Djj,repmat(MD,N,1));
end
