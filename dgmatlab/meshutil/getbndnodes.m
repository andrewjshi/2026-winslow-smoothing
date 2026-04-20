function [e,edg]=getbndnodes(msh,data,bndix)

if nargin<3, bndix=1:-min(msh.t2t(:)); end

ns=size(msh.s,1);
[j,it]=find(ismember(msh.t2t,-bndix));
e=[];
for i=1:length(it)
  newe=ns*(it(i)-1)+double(data.egix(:,j(i)))+1;
  e=[e;newe];
end

evec=zeros(ns*size(msh.t,2),1);
evec(e)=1;
[r,foo,W]=mapdg2cg(msh.p1);
%W(W~=0)=1;
Wevec=W*evec;
e=find(Wevec>0);

evec=W'*Wevec;
edg=find(evec);


