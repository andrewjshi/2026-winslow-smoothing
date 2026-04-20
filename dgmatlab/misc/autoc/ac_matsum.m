function out = ac_matsum(m,n,alpha,A,beta,B,gamma,C)
% Matrix sum
% out = ac_matsum(m,n,alpha,A,beta,B,gamma,C)
%
% C := alpha*A+beta*B+gamma*C
%
% All are m-by-n matrices
if (isnumeric(alpha))
  alpha = num2str(alpha);
end
if (isnumeric(beta))
  beta = num2str(beta);
end
if (isnumeric(gamma))
  gamma = num2str(gamma);
end

out = cell(m,n);
for i=1:m
  for j=1:n
    eqn = [C int2str(i) int2str(j) '='];
    if (~strcmp(alpha,'0'))
      if (~strcmp(alpha,'1'))
        eqn = [eqn alpha '*'];
      end
      eqn = [eqn A int2str(i) int2str(j)];
    end
    if (~strcmp(beta,'0'))
      eqn = [eqn '+'];
      if (~strcmp(beta,'1'))
        eqn = [eqn beta '*'];
      end
      eqn = [eqn B int2str(i) int2str(j)];
    end
    if (~strcmp(gamma,'0'))
      eqn = [eqn '+'];
      if (~strcmp(gamma,'1'))
        eqn = [eqn gamma '*'];
      end
      eqn = [eqn C int2str(i) int2str(j)];
    end
    out{i,j} = eqn;
  end
end
out = reshape(out,m*n,1);
