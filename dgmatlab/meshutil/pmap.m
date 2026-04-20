function p=pmap(p,x,y)
%PMAP  Bilinear mapping on nodes in p (currently only 2-D)
%
%    x,y gives new coordinates of the four unit square corners

s=[0,1,0,1];
t=[0,0,1,1];

C=x/[1,1,1,1; s; t; s.*t];
D=y/[1,1,1,1; s; t; s.*t];

px=p(:,1); py=p(:,2);
p=[0*px+1, px, py, px.*py]*[C;D]';
