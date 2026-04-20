function [q,W,W1,iq]=mapdg2cg(p1,tol)

if nargin<2, tol=1e-8; end

[ns,dim,nt]=size(p1);
p11=reshape(permute(p1,[1,3,2]),[ns*nt,dim]);

snap=max(max(p11,[],1)-min(p11,[],1),[],2)*tol;
[foo,iq,q]=unique(snap*round(p11/snap),'rows');

W1=sparse(q,1:ns*nt,1);
W=spdiags(1./sum(W1,2),0,size(W1,1),size(W1,1))*W1;
