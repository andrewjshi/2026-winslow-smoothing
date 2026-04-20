function msh=tets2msh(XT)
%TETS2MSH  Create msh data structure from block tetrahedra (without bc's)
%
%   Syntax: msh=tets2msh(XT)
%
%   EXAMPLES
%
%   Uniform cube of (5,5,5):
%     dgmeshplot(tets2msh(block2tets(block([5,5,5]))));
%
%   Graded node distribution:
%     dgmeshplot(tets2msh(block2tets(block([20,20,20], -15*ones(12,1)))));
%
%   Boundary layer in 2D:
%     dgmeshplot(tets2msh(block2tets(block([20,20], [1,1,20,20],[0,0;3,0;0,1;3,1]))));

dim = size(XT,1);

p0=reshape(XT,dim,[])';
t0=reshape(1:prod(size(XT))/dim,dim+1,size(XT,3))';
[p,t]=fixmesh(p0,t0);
msh=ml2msh(p,t);
