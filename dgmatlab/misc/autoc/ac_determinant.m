function out = ac_determinant(A, d)
% Determinant
% out = ac_determinant(A,d)
%
% A is 3-by-3

detstr = 'a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31';
out = [d, ' = ', strrep(detstr, 'a', A)];
