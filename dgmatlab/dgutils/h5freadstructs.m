function varargout=h5freadstructs(fname,varargin)

for i=1:nargin-1
  fld=varargin{i};
  fn=sprintf('%s%s.h5',fname,fld);
  varargout{i}=h5freadstruct(fn,dgfieldnames(fld));
end
