function [lM0R,lM0]=mkcholM0(porder,s)

[gx,gw]=gaussquad(2*porder,size(s,2)-1);
gx=[gx,1-sum(gx,2)];
gfs=mkvolshapes(porder,s,gx);

ng=size(gw,1);
lM0=permute(gfs(:,1,:),[1,3,2])*diag(gw,0)* ...
    permute(gfs(:,1,:),[1,3,2])';
lM0R=chol(lM0);
