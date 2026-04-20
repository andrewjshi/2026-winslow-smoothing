function out = ac_array(F,m,n)

if (nargin == 2)
  out = cell(m,1);
  for i=1:m
    out{i} = [F int2str(i)];
  end
elseif (nargin == 3)
  out = cell(m,n);
  for i=1:m
    for j=1:n
      out{i,j} = [F int2str(i) int2str(j)];
    end
  end
  out = reshape(out,m*n,1);
end
