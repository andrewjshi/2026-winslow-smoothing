function out = ac_inner(n,alpha,A,B,beta,C)
% Inner product
% out = ac_inner(n,alpha,A,B,beta,C)
%
% C := alpha*A'*B+beta*C
%
% A is n-by-1, B is n-by-1, C is 1-by-1.

if (isnumeric(alpha))
  alpha = num2str(alpha);
end
if (isnumeric(beta))
  beta = num2str(beta);
end

eqn = [C '='];
if (~strcmp(alpha,'1'))
  eqn = [eqn alpha '*'];
end
eqn = [eqn '('];
for i=1:n
  eqn = [eqn A int2str(i) '*' B int2str(i)];
  if (i ~= n)
    eqn = [eqn '+'];
  end
end
eqn = [eqn ')'];
if (~strcmp(beta,'0'))
  eqn = [eqn '+'];
  if (~strcmp(beta,'1'))
    eqn = [eqn beta '*'];
  end
  eqn = [eqn C];
end

out{1} = eqn;
