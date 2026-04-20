function u1=dginterp(u0,p0,p1,dim)
%DGINTERP  Interpolate, change polynomial order p0->p1
%
%    Syntax: u1=dginterp(u0,p0,p1,[dim]);
%
%    u0:  solution at p0 [ns0,ncomp,nt]
%    p0:  polynomial order p0 (original)
%    p1:  polynomial order p1 (desired)
%    dim: space dimensions (default=2)
%    u1:  solution at p1 [ns1,ncomp,nt]

warning('Cannot find mex-file -- using old ML code (simplex only).');

if nargin<4, dim=2; end

s0=lagrangepnts(p0,dim);
s1=lagrangepnts(p1,dim);

AI=pmonomial(s1(:,1:dim),p0)/pmonomial(s0(:,1:dim),p0);

[ns,ncomp,nt]=size(u0);
u1=reshape(AI*reshape(u0,ns,[]),[],ncomp,nt);
