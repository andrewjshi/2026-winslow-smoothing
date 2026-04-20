function x=lobattopoints(n)

if n==1, x=0; return; end

P=zeros(n+1,n+1);
P([1,2],1)=1;
for k=1:n-1
  P(k+2,1:k+2)=((2*k+1)*[P(k+1,1:k+1) 0]-k*[0 0 P(k,1:k)])/(k+1);
end
x=sort(roots(polyder(P(n,1:n))));
x=[-1;x;1];
