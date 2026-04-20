function msh=mshdiamond1

% Parameters
n1=5;
n2=5;
n3=5;
n4=10;
alpha=2;

% Meshes in parameter space
[p01,e01,t01]=squaremesh(n1,n2);
[p02,e02,t02]=squaremesh(n3,n2);
[p03,e03,t03]=squaremesh(n4,n2);
np01=size(p01,1); nt01=size(t01,1);
np02=size(p02,1); nt02=size(t02,1);
np03=size(p03,1); nt02=size(t03,1);

% Boundary layer
%p01(:,2)=(exp(2*p01(:,2))-1)/(exp(2*1)-1);

% Map
p1=pmap(p01,[0,1,0,1],[.5,0,1.5,1.5]);
p2=pmap(p02,[1,2,1,2],[0,0,1.5,1.5]);
p3=pmap(p03,[2,5,2,5],[0,0,1.5,1.5]);

% Join
p=[p1;p2];
t=[t01;t02+np01];

px=p(:,1); py=p(:,2);
p=[px,py; -px,py; px,-py; -px,-py];
t=[t; t+(np01+np02); t+2*(np01+np02); t+3*(np01+np02)];

t=[t; t03+size(p,1); t03+size(p,1)+np03];
p=[p;p3;p3(:,1),-p3(:,2)];

% Remove duplicated nodes and orient elements
[p,t]=fixmesh(p,t);

% Create boundary numbering and msh structure
bndexpr={'all(p(:,1).^2/1.1^2+p(:,2).^2/.7^2<1)','any(p(:,1).^2/1.1^2+p(:,2).^2/.7^2>=1)'};
msh=ml2msh(p,t,bndexpr);
msh=mshcurved(msh,[]);
