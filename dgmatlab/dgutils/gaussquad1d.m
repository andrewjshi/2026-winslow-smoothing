function [x,w]=gaussquad1d(p)
%GAUSSQUAD1D  Gaussian quadrature on [0,1] for given degree of precision

n = ceil((p+1)/2);
b = 1:n-1;
b = b ./ sqrt(4*b.^2-1);
[Q,D] = eig(diag(b,1) + diag(b,-1));
x = (diag(D) + 1) / 2;
w = (Q(1,:).^2)';
