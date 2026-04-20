function D=blockdiag(A,n)

[i,j,s]=find(A);
blk=floor((i-1)/n);
i1=n*blk+1;
i2=n*blk+n;
ix=j>=i1 & j<=i2;
D=sparse(i(ix),j(ix),s(ix),size(A,1),size(A,2));
