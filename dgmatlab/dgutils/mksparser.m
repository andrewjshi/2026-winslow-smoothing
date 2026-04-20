function A=mksparser(A)
% Fix up sparse matrix if e.g. Matlab allocates too much memory

[i,j,s]=find(A);
A=sparse(i,j,s,size(A,1),size(A,2));
