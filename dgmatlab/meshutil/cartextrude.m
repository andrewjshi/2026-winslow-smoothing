function [x,y,z]=cartextrude(x,y,nz,Lz)

% Extrude in z-direction
z=repmat(permute(linspace(0,Lz,nz),[1,3,2]),[size(x),1]);
x=repmat(x,[1,1,nz]);
y=repmat(y,[1,1,nz]);

% Flip y/z to make x/z the symmetry plane
foo=y; y=z; z=foo;
x=permute(x,[1,3,2]);
y=permute(y,[1,3,2]);
z=permute(z,[1,3,2]);
