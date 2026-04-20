function [p,t]=block2pq(X)

XT=block2hexes(X);
dim = size(XT,1);
p0=reshape(XT,dim,[])';
t0=reshape(1:prod(size(XT))/dim,2^dim,size(XT,3))';
[p,t]=fixmesh(p0,t0);
