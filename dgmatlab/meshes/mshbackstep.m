function msh=mshbackstep

X=block([20,20],[10,1,-10,1/10],[0,0;1,0;0,1;1,1]);
Y=block([20,20],[10,10,-10,10],[0,-1;1,-1;0,0;1,0]);

[p1,t1]=block2pt(X);
[p2,t2]=block2pt(Y);

p=[p1;p2];
t=[t1;t2+size(p1,1)];

[p,t]=fixmesh(p,t);

[t2t,t2n]=mkt2t(t);
msh=ml2msh(p,t,t2t,t2n);
msh=mshcurved(msh,[]);
