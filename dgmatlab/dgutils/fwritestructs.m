function fwritestructs(fname,varargin)

for i=2:nargin
  fld=inputname(i);
  fn=sprintf('%s%s.dat',fname,fld);
  fwritestruct(varargin{i-1},fn,dgfieldnames(fld));
end
