function Dinv=invblockdiag(A,k)

n=size(A,1);
nt=n/k;
[ii,jj,kk]=ndgrid(1:k,1:k,1:nt);
iii=ii+(kk-1)*k;
jjj=jj+(kk-1)*k;

sss=zeros(k,k,nt);
for it=1:nt
  i=k*(it-1)+(1:k);
  sss(:,:,it)=inv(full(A(i,i)));
end

Dinv=sparse(iii(:),jjj(:),sss(:),n,n);
