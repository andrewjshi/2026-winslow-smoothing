function msh=ml2msh(p,t,bndexpr,fd,fdargs,periodic_map,curvedbnds)

if nargin<3, bndexpr=[]; end
if nargin<4, fd=[]; end
if nargin<5, fdargs={}; end
if nargin<6, periodic_map=[]; end
if nargin<7, curvedbnds=[]; end

eltype=dimnv2eltype(size(p,2),size(t,2));
[t2t,t2n]=mkt2t(t,true,eltype);
if ~isempty(bndexpr)
  t2t=setbndnbrs(p,t,t2t,bndexpr);
end
if ~isempty(periodic_map)
  [t2t,t2n]=setbndperiodic(p,t,t2t,t2n,periodic_map);
end

msh.p=p';
msh.t=int32(t'-1);
msh.t2t=int32(t2t'-1);
msh.t2n=int32(t2n'-1);
msh.bndexpr=bndexpr;
msh.fd=fd;
msh.fdargs=fdargs;
msh.eltype=int32(eltype);

msh=mshcurved(msh,curvedbnds);
