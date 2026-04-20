function out = gemv(trans, m, n, alpha, A, X, beta, Y)
% GEMV  performs one of the matrix-vector operations
% out = gemv(trans, m, n, alpha, A, X, beta, Y)
%
%    y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
% 
% where alpha and beta are scalars, x and y are vectors and A is an
% m by n matrix.
%
% trans is 'n' to use A, and 't' to use A'.
%
% For example:
% ac_gemv('n', 3, 3, 1, 'G', 'N', 0, 'n')
%  'n1=(G11*N1+G12*N2+G13*N3)'
%  'n2=(G21*N1+G22*N2+G23*N3)'
%  'n3=(G31*N1+G32*N2+G33*N3)'

if (trans == 'N' || trans == 'n')
  t = false;
elseif (trans == 'T' || trans == 't')
  % Swap m and n
  t = n;
  n = m;
  m = t;
  t = true;
else
  warning('invalid trans: allowed values are n and t');
end

if (isnumeric(alpha))
  alpha = num2str(alpha);
end
if (isnumeric(beta))
  beta = num2str(beta);
end

out = cell(m,1);
for i=1:m
  eqn = [Y, int2str(i), '=('];
  for j=1:n
    eqn = [eqn A];
    if (t)
      eqn = [eqn int2str(j) int2str(i)];
    else
      eqn = [eqn int2str(i) int2str(j)];
    end
    eqn = [eqn '*' X int2str(j)];
    if (j ~= n)
      eqn = [eqn '+'];
    end
  end
  eqn = [eqn ')'];
  if (~strcmp(alpha,'1'))
    eqn = [eqn '*' alpha];
  end
  if (~strcmp(beta,'0'))
    eqn = [eqn '+' beta '*' Y int2str(i)];
  end
  out{i,1} = eqn;
end
