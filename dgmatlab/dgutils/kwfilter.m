function ulow=kwfilter(msh,plow,u)

p = double(msh.porder);
d = size(msh.p,1);

nnodes = [2 3 4 5 6;3 6 10 15 21;4 10 20 35 56];

nlow = nnodes(d,plow);
nhigh = nnodes(d,p) - nnodes(d,plow);

V = kwvander(p,2*msh.s(:,1:d)-1);
PL = diag([ones(nlow,1); zeros(nhigh,1)]);
F = V*PL*inv(V);

ulow = reshape(F*reshape(u,size(u,1),[]),size(u));
