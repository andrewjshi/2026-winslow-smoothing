function t2t=setbndnbrs(p0,t,t2t,bndexpr)

dim=size(p0,2);
nf=size(t2t,2);
eltype=dimnv2eltype(size(p0,2),size(t,2));
facemap=mkfacemap(dim,eltype);

[i,j]=find(t2t==0);
for ii=1:length(i)
  tix=facemap(:,j(ii));
  n=t(i(ii),tix);
  p=p0(n,:);

  found=false;
  for jj=1:length(bndexpr)
    if ischar(bndexpr{jj})
        bndval = eval(bndexpr{jj});
    else
        bndval = bndexpr{jj}(p);
    end
    if bndval
      found=true;
      bnd=jj;
      break;
    end
  end
  
  if ~found
    error('Strange boundary.');
  end
  
  t2t(i(ii),j(ii))=-bnd+1;
end
