function viz = dgvtk(fname, msh, u, navstok)
%DGVTL Export DG data to VTK / Paraview
%
%   dgvtk(fname, msh, u=[], navstok=false)
%
%   Only provided for backward compatibility - use DGVTKWRITE instead.

    if nargin < 2, u = []; end
    if nargin < 3, navstok = 0; end
    dgvtkwrite(fname, dgviz_vol(msh, u, navstok));
    
end
