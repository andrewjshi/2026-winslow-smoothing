function e=getcornernodes(msh)

ns=size(msh.s,1);
nt=size(msh.t,2);

ix=find(sum(msh.s<1e-3,2)==2);
e=ns*(repmat(1:nt,3,1)-1)+repmat(ix(:),1,nt);

evec=zeros(ns*size(msh.t,2),1);
evec(e)=1;
[r,W]=mapdg2cg(msh.p1);
W(W~=0)=1;
Wevec=W*evec;
e=find(Wevec>0);
