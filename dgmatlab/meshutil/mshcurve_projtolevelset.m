function msh=mshcurve_projtolevelset(msh,fd,fdargs)

for it=1:size(msh.t,2)
  p1=msh.p1(:,:,it);
  p=msh.p(:,msh.t(:,it)+1)';
  d=feval(fd,p,fdargs{:});
  d0=msh.s*d;
  for iii=1:5
    p1=pproj(p1,d0,fd,fdargs{:});
  end
  msh.p1(:,:,it)=p1;
end
msh = mshcurved(msh, 'all');
