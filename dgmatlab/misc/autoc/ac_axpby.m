function out = ac_axpby(n, alpha, X, beta, Y)
% AXPY  constant times a vector plus a constant times a vector.
% out = axpby(n, alpha, X, beta, Y)
%
%    y := alpha*x + beta*y
% 
% where alpha and beta are scalars, x and y are vectors of length n
%
% For example:
% ac_axpby(2,'alpha','X','beta','Y')
%  'Y1=alpha*X1+beta*Y1'
%  'Y2=alpha*X2+beta*Y2'

if (isnumeric(alpha))
  alpha = num2str(alpha);
end
if (isnumeric(beta))
  beta = num2str(beta);
end

out = cell(n,1);
for i=1:n
  eqn = [Y int2str(i) '='];
  if (strcmp(alpha,'0'))
    eqn = [eqn '0'];
  elseif (strcmp(alpha,'1'))
    eqn = [eqn X int2str(i)];
  else
    eqn = [eqn X int2str(i) '*(' alpha ')'];
  end
  if (~strcmp(beta,'0'))
    eqn = [eqn '+'];
    eqn = [eqn Y int2str(i)];
    if (~strcmp(beta,'1'))
      eqn = [eqn '*(' beta ')'];
    end
  end
  out{i} = eqn;
end
