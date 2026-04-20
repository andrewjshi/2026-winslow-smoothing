function out = ac_scal(n,alpha,A,B)
% SCAL
% out = ac_scal(n,alpha,a,[b])
%
% b := alpha*a
% where a and b are vectors of length n
%
% if b is not provided, do
% a := alpha*a

if (isnumeric(alpha))
  alpha = num2str(alpha);
end

if (nargin<4)
  B = A;
end

out=cell(n,1);

for i=1:n
  out{i} = [B int2str(i) '=(' alpha ')*' A int2str(i)];
end

  