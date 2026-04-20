function [ii,jj] = dgindices_mass(ns,nt)

j = repmat(1:ns,ns,1);
i = j';
jj = repmat(j(:),1,nt) + ns*(0:nt-1);  
ii = repmat(i(:),1,nt) + ns*(0:nt-1);  
