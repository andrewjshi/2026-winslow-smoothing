function v=dgvorticity(msh,u)

switch size(msh.p,1)
 case 2
  uv=u(:,2:3,:)./u(:,[1,1],:);
  ux=dgdiff(msh,uv,1);
  vx=dgdiff(msh,uv,2);
  
  v=ux(:,2,:)-vx(:,1,:);
 case 3
  uvw=u(:,2:4,:)./u(:,[1,1,1],:);
  ux=dgdiff(msh,uvw,1);
  vx=dgdiff(msh,uvw,2);
  wx=dgdiff(msh,uvw,3);
  
  v=sqrt((vx(:,3,:)-wx(:,2,:)).^2 + ...
         (ux(:,3,:)-wx(:,1,:)).^2 + ...
         (ux(:,2,:)-vx(:,1,:)).^2);
end
