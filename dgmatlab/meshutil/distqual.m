function ga=distratio(msh,bnds,nneighbors,chunk)

if nargin<2, bnds=1:-min(msh.t2t(:)); end
if nargin<3, nneighbors=0; end
if nargin<4, chunk=inf; end

if isempty(nneighbors) | nneighbors<0
  elix=1:size(msh.t,2);
else
  elix=mkbndneighbors(msh,bnds,nneighbors);
end

if isinf(chunk)
  J=geojac(msh);
else
  ct=1;
  nt=size(msh.p1,3);
  J=[];
  while ct<=nt
    lt=min(ct+chunk-1,nt);
    ix=ct:lt;
    msh1=msh;
    msh1.p1=msh.p1(:,:,ix);
    J=[J,geojac(msh1)];
    ct=lt+1;
  end
end

minJ=min(J,[],1);
maxJ=max(J,[],1);

ga=minJ(elix)./maxJ(elix);
