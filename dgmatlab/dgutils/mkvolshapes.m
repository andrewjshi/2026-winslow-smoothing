function gfs=mkvolshapes(porder,s,gx)
%MKVOLSHAPES  Shape functions and derivatives at Gauss points
%
%    Syntax: gfs=mkvolshapes(porder,s,gx);
%
%    porder: polynomial order
%    s:      node points in element (barycentric)
%            dimensions [npoints,dim+1]
%    gx:     gauss points in element (barycentric)
%            dimensions [ngauss,dim+1]
%    gfs:    shape functions and derivatives
%            dimensions [npoints,dim+1,ngauss]
%            gfs(:,1,:) values
%            gfs(:,1+d,:cat ') derivatives wrt coordinate d

dim=size(s,2)-1;
[gf,gfx,gfy,gfz]=pmonomial(gx(:,2:dim+1),porder);
A=pmonomial(s(:,2:dim+1),porder);
gfs=gf/A;
if dim>=1, gfs=[gfs,gfx/A]; end
if dim>=2, gfs=[gfs,gfy/A]; end
if dim>=3, gfs=[gfs,gfz/A]; end
if dim>=4, error('Dimension not implemented.'); end
gfs=reshape(gfs,[size(gf),dim+1]);
gfs=permute(gfs,[2,3,1]);
