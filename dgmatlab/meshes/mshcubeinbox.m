function msh=mshcubeinbox(mul1)

if nargin<1, mul1=1; end

m=3*mul1+1;
n=3*2*mul1+1;

[p,e,t]=cubemesh(n,m,m);

p(:,[2,3])=p(:,[2,3])*6-3;
p(:,1)=p(:,1)*12-3;

pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:)+p(t(:,4),:))/4;
t(all(abs(pmid)<1-1e-6,2),:)=[];

[p,t]=fixmesh(p,t);

bndexpr={'all(sum(p.^2,2)<2^2)','any(sum(p.^2,2)>=2^2)'};
msh=ml2msh(p,t,bndexpr);
msh=mshcurved(msh,[]);
