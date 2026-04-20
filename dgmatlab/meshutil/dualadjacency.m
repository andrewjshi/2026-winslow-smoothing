function A=dualadjacency(t2t,val)

nf=size(t2t,2);
nt=size(t2t,1);

row=(1:nt)'*ones(1,nf);
col=double(t2t);
col(t2t<=0)=1;
if nargin<2
  val=ones(nt,nf);
end
val(t2t<=0)=0;

A=sparse(row,col,val,nt,nt);
