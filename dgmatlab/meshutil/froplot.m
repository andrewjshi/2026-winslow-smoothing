function froplot(fro,R0)

if nargin<1, R0=20; end

p=fro.p;
e=fro.e(:,1:3);
en=fro.e(:,4);

pmid=(p(e(:,1),:)+p(e(:,2),:)+p(e(:,3),:))/3;
rmid=sqrt(sum(pmid.^2,2));
out=unique(en(rmid>R0));
keep=setdiff(1:max(en),out);
ekeep=ismember(en,keep);

clf
for ikeep=keep
  ce=e(en==ikeep,:);
  n=0*p;
  [foo,cn]=froseval(fro,ikeep,fro.ps{ikeep}(:,2:3));
  n(fro.ps{ikeep}(:,1),:)=cn;
  hh=patch('vertices',p,'faces',ce,'vertexnormals',n);
  set(hh,'edgecolor','k','facecolor',[.8,1,.8]);
end

view(3),axis equal,axis off
lighting gouraud,camlight(0,0)
fancycamera
