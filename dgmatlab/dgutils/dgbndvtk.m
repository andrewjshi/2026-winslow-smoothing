function dgbndvtk(fname, msh, u, app, phys, postqtys, bndnbrs)
%DGBNDVTK Export boundary DG data to VTK / Paraview
%
%   dgbndvtk(fname, msh, u, app, phys, postqty, bndnbrs='all')
%
%   Only provided for backward compatibility - use DGVTKWRITE instead.

    if nargin < 7, bndnbrs = 'all'; end

    dgvtkwrite(fname, dgviz_bnd(msh, u, app, phys, postqtys, bndnbrs));
    
end
