function h=dgplot(varargin)
%DGPLOT  Plot DG solution
%
%    Syntax: dgplot(msh,u,[utype,clim,nref,pltmesh,bnds])
%
%    msh:     mesh structure
%    u:       solution
%    utype:   physical quantity to plot (see dgeval,default=component 1)
%    clim:    limits on color or y scale (see MATLAB manual, default=[])
%    nref:    number of uniform refinements of original mesh for plotting
%             (default=ceil(log2(porder)))
%    pltmesh: plot mesh or not (logical or color, default=false)
%    bnds:    boundary surfaces to plot (3-D only, default all boundaries)

dim=size(varargin{1}.p,1);
switch dim
 case 1
  h=dgplot1d(varargin{:});
 case 2
  dgplot2d(varargin{:});
 case 3
  h=dgplot3d(varargin{:});
 otherwise('Unknown dimension.');
end

if nargout<1, clear h; end
