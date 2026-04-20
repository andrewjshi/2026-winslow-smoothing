function msh=mk3Dcartfoilmesh(porder,xm,ym,zm,bndexpr)

mid=(size(xm,3)+1)/2;
[pa,ta,p1a]=cart2dg(porder,xm(:,:,1:mid),ym(:,:,1:mid),zm(:,:,1:mid),struct('ind',0));
[pb,tb,p1b]=cart2dg(porder,xm(:,:,mid:end),ym(:,:,mid:end),zm(:,:,mid:end),struct('ind',2));
p=[pa;pb]; t=[ta;tb+size(pa,1)]; p1=cat(3,p1a,p1b);
[p,t]=fixmesh(p,t,1e-8);
msh=cart2msh(porder,p,t,p1,bndexpr);
msh=cartmsh_yperiodic(msh);
