function [y,dy]=plegendre(x,n1)

y=zeros(length(x),n1+1);
dy=zeros(length(x),n1+1);
y(:,1)=1;
dy(:,1)=0;
if n1>=1
  y(:,2)=x;
  dy(:,2)=1;
end
for ii=1:n1-1
  y(:,ii+2)=((2*ii+1)*x.*y(:,ii+1)-ii*y(:,ii))/(ii+1);
  dy(:,ii+2)=((2*ii+1)*x.*dy(:,ii+1)+(2*ii+1)*y(:,ii+1)-ii*dy(:,ii))/(ii+1);
end
