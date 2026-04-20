function msh=mshdiamond(alpha,nscale)

% Parameters
if nargin<1, alpha=4; end
if nargin<2, nscale=2; end

n1=5*nscale;
n2=10*nscale;
n3=6*nscale;
n4=2*nscale;

% Meshes in parameter space
[p01,e01,t01]=squaremesh(n1,n2);
[p02,e02,t02]=squaremesh(n3,n1);
[p03,e03,t03]=squaremesh(n4,2*n1-1);
np01=size(p01,1); nt01=size(t01,1);
np02=size(p02,1); nt02=size(t02,1);
np03=size(p03,1); nt03=size(t03,1);

% Boundary layer
p01(:,2)=(exp(alpha*p01(:,2))-1)/(exp(alpha*1)-1);

% Map
p1=pmap(p01,[0,1,0,5],[.5,0,2.5,0]);
p2=pmap(p02,[5,7,0,7],[0,0,2.5,2.5]);
p3=pmap(p03,[7,8,7,8],[-2.5,-2.5,2.5,2.5]);

% Join
p=[p1;p2];
t=[t01;t02+np01];

px=p(:,1); py=p(:,2);
p=[px,py; -px,py; px,-py; -px,-py];
t=[t; t+(np01+np02); t+2*(np01+np02); t+3*(np01+np02)];

t=[t; t03+size(p,1)];
p=[p;p3];

% Remove duplicated nodes and orient elements
[p,t]=fixmesh(p,t);

% Create boundary numbering and msh structure
bndexpr={'all(p(:,2)<-2.5+1e-6)','all(p(:,1)>8-1e-6)', ...
         'all(p(:,2)>2.5-1e-6)','all(p(:,1)<-7+1e-6)', ...
         'all(sum(p.^2,2)<2^2)'};
msh=ml2msh(p,t,bndexpr);
msh=mshcurved(msh,[]);
