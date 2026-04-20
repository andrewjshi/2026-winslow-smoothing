function varargout=freadstructs(fname,varargin)

for i=1:nargin-1
  fld=varargin{i};
  fn=sprintf('%s%s.dat',fname,fld);
  varargout{i}=freadstruct(fn,dgfieldnames(fld));
end
