function msh=dgnodesfro(msh,nsmooth)

if nargin<2, nsmooth=0; end

porder=msh.porder;
dim=size(msh.p,1);

[q,W,W1]=mapdg2cg(msh.p1,1e-5);
q=reshape(q,size(msh.p1,1),size(msh.p1,3));

cp1=W*reshape(permute(msh.p1,[1,3,2]),[],3);

[is,js]=find(msh.fromap);
for ix=1:length(is)
  i=is(ix); j=js(ix);
  six=find(msh.s(:,j)<1e-6);
  spnts=setdiff(1:4,j);
  
  cie=msh.fromap(i,j);
  [foo,map]=ismember(msh.t(spnts,i)+1,msh.fro.e(cie,1:3));
  [foo,map]=sort(map);
  cs=msh.s(six,spnts(map));
  
  pnew=frotriinterp(msh.fro,cie,cs);
  for iter=1:nsmooth
    csnew=sequiv(pnew,cs);
    pnew=frotriinterp(msh.fro,cie,csnew);
  end

  cp1(q(six,i),:)=pnew;
end

newp1=permute(reshape(W1'*cp1,size(msh.p1,1),size(msh.p1,3),3),[1,3,2]);

% Todo: Figure out real curved
msh.ncurved(:)=true;
msh.tcurved(:)=true;
msh.ecurved(:)=true;

msh.p1=newp1;
