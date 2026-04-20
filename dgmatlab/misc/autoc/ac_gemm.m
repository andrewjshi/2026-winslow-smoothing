function out = ac_gemv(transA, transB, m, n, k, alpha, A, B, beta, C)
% DGEMM  performs one of the matrix-matrix operations
% out = ac_gemv(transA, transB, m, n, k, alpha, A, B, beta, C)
%
%    C := alpha*op( A )*op( B ) + beta*C,
% 
% where  op( X ) is one of
% 
%    op( X ) = X   or   op( X ) = X',
% 
% alpha and beta are scalars, and A, B and C are matrices, with op( A )
% an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
%
% ac_gemm('n', 'n', 3,1,2,1,'A','B',0,'C')
%   'C11=(A11*B11+A12*B21)'
%   'C21=(A21*B11+A22*B21)'
%   'C31=(A31*B11+A32*B21)'

if (transA == 'N' || transA == 'n')
  tA = false;
elseif (transA == 'T' || transA == 't')
  tA = true;
else
  warning('invalid transA: allowed values are n and t');
end

if (transB == 'N' || transB == 'n')
  tB = false;
elseif (transB == 'T' || transB == 't')
  tB = true;
else
  warning('invalid transB: allowed values are n and t');
end

if (isnumeric(alpha))
  alpha = num2str(alpha);
end
if (isnumeric(beta))
  beta = num2str(beta);
end

out = cell(m,n);
for i=1:m
  for j=1:n
    eqn = [C int2str(i) int2str(j) '='];
    if(~strcmp(alpha,'1'))
      eqn = [eqn alpha '*'];
    end
    eqn = [eqn '('];
    for s=1:k
      if (tA)
        eqn = [eqn A int2str(s) int2str(i)];
      else
        eqn = [eqn A int2str(i) int2str(s)];
      end
      eqn = [eqn '*'];
      if (tB)
        eqn = [eqn B int2str(j) int2str(s)];
      else
        eqn = [eqn B int2str(s) int2str(j)];
      end
      if (s~=k)
        eqn = [eqn '+'];
      end
    end
    eqn = [eqn ')'];
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

out = reshape(out,i*j,1);
