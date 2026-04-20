function d=mshdist(msh,bnds,p1)
%MSHDIST  Compute distance function for a mesh.
%   d=mshdist(msh,bnds)
%     msh=msh from nodealloc
%     bnds=[list of boundary numbers]
%     d=smallest absolute distance from node to bnds
%
%   d=mshdist(msh,bnds,p1)
%     p1=nodes to evaluate distance at (instead of msh.p1)
