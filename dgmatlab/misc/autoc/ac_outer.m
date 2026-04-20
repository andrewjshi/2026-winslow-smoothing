function out = ac_outer(m,n,alpha,A,B,beta,C,D)
% Outer product
% out = ac_outer(m,n,alpha,A,B,beta,C,D)
%
% D := alpha*A*B'+beta*C
%
% A is m-by-1, B is n-by-1, C is m-by-n, D is m-by-n.
%
% D is optional, and will be replaced by C if it is not provided.

if (isnumeric(alpha))
  alpha = num2str(alpha);
end
if (isnumeric(beta))
  beta = num2str(beta);
end

if (nargin<8)
  D = C;
end

out = cell(m,n);
for i=1:m
  for j=1:n
    eqn = [D int2str(i) int2str(j) '='];
    if (~strcmp(alpha,'1'))
      eqn = [eqn alpha '*'];
    end
    eqn = [eqn A int2str(i) '*' B int2str(j)];
    if (~strcmp(beta,'0'))
      eqn = [eqn '+'];
      if (~strcmp(beta,'1'))
        eqn = [eqn beta '*'];
      end
      eqn = [eqn C int2str(i) int2str(j)];
    end
    out{i,j} = eqn;
  end
end
out = reshape(out,m*n,1);
