function h5fwritestructs(fname,varargin)

for i=2:nargin
  fld=inputname(i);
  fn=sprintf('%s%s.h5',fname,fld);
  h5fwritestruct(varargin{i-1},fn,dgfieldnames(fld));
end
