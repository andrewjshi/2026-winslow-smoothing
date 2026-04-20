function out = ac_norm2(n,A,b)
% 2-Norm
% out = ac_outer(n,A,b)
%
% b := ||A_2|| = sqrt(A1^2+A2^2+...)
%
% A is n-by-1.

eqn = [b '=sqrt('];
for i=1:n
  eqn = [eqn A int2str(i) '*' A int2str(i)];
  if (i ~= n)
    eqn = [eqn '+'];
  end
end
eqn = [eqn ')'];

out{1} = eqn;
