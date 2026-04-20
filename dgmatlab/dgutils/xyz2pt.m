function [p,t,u] = xyz2pt(xyz, uxyz)
%XYZ2PT  Convert mesh in xyz format (+data) to p,t format (+node data)
%
%   [p,t,u] = xyz2pt(xyz, uxyz)
%
%   Dimensions:
%      xyz:  dim x nv x nelem
%      uxyz: nv x nelem
%      p:    np x dim
%      t:    nelem x nv
%      u:    np
    
    [dim,nv,nel] = size(xyz);
    
    t = reshape(1:nv*nel, nv, nel)';
    p = reshape(xyz, dim, nv*nel)';
    
    [p,t,pix] = fixmesh(p,t);
    u = uxyz(pix);
    
end
