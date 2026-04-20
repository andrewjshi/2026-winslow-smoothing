function p=mapp(p,x)
%MAPP  ndim-linear map points
%
%    x(1:2^ndim,1:ndim) gives new coordinates of the unit cube

switch size(p,2)
 case 1
     r=[0,1];
     
     C=x(:,1)'/[1,1; r];
     
     p=[0*p+1,p]*C';
 case 2    
     r=[0,1,0,1];
     s=[0,0,1,1];
     
     C=x(:,1)'/[1,1,1,1; r; s; r.*s];
     D=x(:,2)'/[1,1,1,1; r; s; r.*s];

     px=p(:,1); py=p(:,2);
     p=[0*px+1, px, py, px.*py]*[C;D]';
 case 3
     r=[0,1,0,1,0,1,0,1];
     s=[0,0,1,1,0,0,1,1];
     t=[0,0,0,0,1,1,1,1];
     
     C=x(:,1)'/[1,1,1,1,1,1,1,1; r; s; t; r.*s; r.*t; s.*t; r.*s.*t];
     D=x(:,2)'/[1,1,1,1,1,1,1,1; r; s; t; r.*s; r.*t; s.*t; r.*s.*t];
     E=x(:,3)'/[1,1,1,1,1,1,1,1; r; s; t; r.*s; r.*t; s.*t; r.*s.*t];

     px=p(:,1); py=p(:,2); pz=p(:,3);
     p=[0*px+1, px, py, pz, px.*py, px.*pz, py.*pz, px.*py.*pz]*[C;D;E]';
end
