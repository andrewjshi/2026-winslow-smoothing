function [p,q]=quniref(p,q,nref)

if nargin<3, nref=1; end

for iref=1:nref
  np=size(p,1);
  nq=size(q,1);
  q=q(:,[1,2,4,3]);

  pair=[q(:,[1,2]);q(:,[2,3]);q(:,[3,4]);q(:,[4,1])];
  [pair,pairi,pairj]=unique(sort(pair,2),'rows');
  pmid=(p(pair(:,1),:)+p(pair(:,2),:))/2;
  pc=(p(q(:,1),:)+p(q(:,2),:)+p(q(:,3),:)+p(q(:,4),:))/4;
  q1=q(:,1);
  q2=q(:,2);
  q3=q(:,3);
  q4=q(:,4);
  q12=pairj(1:nq)+np;
  q23=pairj(nq+1:2*nq)+np;
  q34=pairj(2*nq+1:3*nq)+np;
  q41=pairj(3*nq+1:4*nq)+np;
  qc=(1:nq)'+np+size(pmid,1);
  
  q=[q1,q12,qc,q41;
     q2,q23,qc,q12;
     q3,q34,qc,q23;
     q4,q41,qc,q34];
  p=[p;pmid;pc];
  q=q(:,[1,2,4,3]);
end
