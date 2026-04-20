function M = dgmass_ml(msh, data)

phi = permute(data.gfs(:,1,:),[1,3,2]);
J = geojac(msh, data);

[ns,dim,nt] = size(msh.p1);

MM=zeros(ns,ns,nt);
for it = 1:nt
  mul = data.gw.*J(:,it)/2;
  M = phi*diag(mul)*phi';
  MM(:,:,it) = M;
end

[ii,jj] = dgindices_mass(ns, nt);
[r,W] = mapdg2cg(msh.p1);
M = sparse(r(ii(:)),r(jj(:)),MM(:));
