function msh=mshfix(msh)
%MSHFIX  NOTE: New functionality, order for ILU performance

% Old functionality:
%MSHFIX  Move curved elements first.

msh=mshreorder(msh, @symrcm);
return;

%Old functionality:

if ~isfield(msh,'ecurved')
  error('No ecurved, run mshcurved.');
end
if ~isfield(msh,'tcurved')
  msh.tcurved=any(msh.ecurved,1);
end

perm=[find(msh.tcurved(:)); find(~msh.tcurved(:))];
msh=mshreorder(msh,perm);
