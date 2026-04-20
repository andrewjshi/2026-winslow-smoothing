function a = kwvander(n,xn)

switch size(xn,2)
 case 1
  a = kwvander1d(n,xn);
 case 2
  a = kwvander2d(n,xn);
 case 3
  a = kwvander3d(n,xn);
 otherwise
  error('Dimension not implemented.');
end

function a = kwvander1d(n,xn)

for i=1:(n+1)
  pp = jacobi(i-1,0,0);
  a(:,i) = polyval(pp,xn(:));
end
        
function a = kwvander2d(n,xn)

pq = pascalindex2d(n);
for i=1:(n+1)*(n+2)/2;
  a(:,i) = koornwinder2d(xn,pq(i,:));
end

function a = kwvander3d(n,xn)

pq = pascalindex3d(n);
for i=1:(n+1)*(n+2)*(n+3)/6;
  a(:,i) = koornwinder3d(xn,pq(i,:));
end


function pq = pascalindex3d(n)

l=1;
for i=0:n
  for j=0:i
    for k=0:j
      pq(l,1)=i-j;
      pq(l,2)=j-k;
      pq(l,3)=k;
      l = l+1;
    end
  end
end

function pq = pascalindex2d(n)

l=1;
for i=0:n
  for j=0:i
    pq(l,1)=i-j;
    pq(l,2)=j;
    l = l+1;
  end
end

function x = koornwinder2d(xc,pq)
%x = koornwinder2d(xc,pq)
%Return the values of the Triangular Koornwinder 
%Polynomial pq for the input values xc

for i=1:length(xc(:,1))
  if (xc(i,2) < 1) 
    e(i,1) = 2*(1+xc(i,1))/(1-xc(i,2))-1;
  else
    e(i,1) = 1;
  end
end
e(:,2) = xc(:,2);

pp = jacobi(pq(1),0,0);
qp = jacobi(pq(2),2*pq(1)+1,0);
for i=1:pq(1)
  qp = conv([-0.5,0.5],qp);
end

pval = polyval(pp,e(:,1));
qval = polyval(qp,e(:,2));

x = pval.*qval;

function x = koornwinder3d(xc,pq)
%x = koornwinder3d(xc,pq)
%Return the values of the Tetrahedral Koornwinder 
%Polynomial pq for the input values xc

for i=1:length(xc(:,1))
  if (xc(i,2)+xc(i,3) < 0)
    e(i,1) = -2*(1+xc(i,1))/(xc(i,2)+xc(i,3)) - 1;
  else
    e(i,1) = -1;
  end
  
  if (xc(i,3) < 1) 
    e(i,2) = 2*(1+xc(i,2))/(1-xc(i,3))-1;
  else
    e(i,2) = 1;
  end
end
e(:,3) = xc(:,3);

pp = jacobi(pq(1),0,0);
qp = jacobi(pq(2),2*pq(1)+1,0);
rp = jacobi(pq(3),2*pq(1)+2*pq(2)+2,0);
for i=1:pq(1)
  qp = conv([-0.5,0.5],qp);
end
for i=1:pq(1)+pq(2)
  rp = conv([-0.5,0.5],rp);
end

pval = polyval(pp,e(:,1));
qval = polyval(qp,e(:,2));
rval = polyval(rp,e(:,3));

x = pval.*qval.*rval;

function p = jacobi(n,a,b)
%p = jacobi(n,a,b)
%Return Jacobi Polynomial (a,b) of order n

p0 = [1];
p1 = [(a+b+2)/2, (a-b)/2];

if n==0
  p = p0;
elseif n==1
  p = p1;
elseif n>1        
  for i=1:n-1
    a1 = 2*(i+1)*(i+a+b+1)*(2*i+a+b);
    a2 = (2*i+a+b+1)*(a*a-b*b);
    a3 = (2*i+a+b)*(2*i+a+b+1)*(2*i+a+b+2);
    a4 = 2*(i+a)*(i+b)*(2*i+a+b+2);
    p2 = conv([a3,a2],p1);
    q = zeros(1,i+2);
    q(3:i+2) = a4*p0;
    p = (p2-q)/a1;
    p0 = p1;
    p1 = p;
  end
end
