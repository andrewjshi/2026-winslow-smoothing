function msh=mshcurved(msh,curvedbnds)
%MSHCURVED  Compute ecurved and ncurved

[dim,np]=size(msh.p);
[nf,nt]=size(msh.t2t);
if ischar(curvedbnds) & strmatch(curvedbnds,'all')
    msh.ncurved=logical(ones(1,np));
    msh.ecurved=logical(ones(nf,nt));
    msh.tcurved=logical(ones(1,nt));
    return;
elseif isempty(curvedbnds)
    msh.ncurved=logical(zeros(1,np));
    msh.ecurved=logical(zeros(nf,nt));
    msh.tcurved=logical(zeros(1,nt));
    return;
end

t=double(msh.t'+1);
t2t=double(msh.t2t'+1);

edges=mkfaces(t,msh.eltype);

if dim==1 | dim==2
  ecurved=ismember(-msh.t2t,curvedbnds);
  ix=find(ecurved');
  ncurvedix=unique(edges(ix,:));
  ncurved=logical(zeros(1,size(msh.p,2)));
  ncurved(ncurvedix)=true;
elseif msh.eltype==t_simplex
  bndnodes=unique(surftri(msh.p',t));
  faces=[t(:,[2,3,4]);
         t(:,[3,4,1]);
         t(:,[4,1,2]);
         t(:,[1,2,3])];
  facenbr=repmat(t2t(:),[1,3]);
  ncurved=~~accumarray(faces(:),double(ismember(-facenbr(:)+1,curvedbnds)))';
  isbnd=ismember(faces,bndnodes) & ismember(faces,find(ncurved));
  ecurved=reshape(sum(isbnd,2)>=2,size(t))';
else
    error('not implemented');
end
msh.ncurved=ncurved;
msh.ecurved=ecurved;
msh.tcurved=any(msh.ecurved,1);
