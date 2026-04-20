function gw=gaussweights(gx)
    
n=size(gx,1);
V=plegendre(gx,n);
f=[2;zeros(n,1)];
gw=V'\f;
