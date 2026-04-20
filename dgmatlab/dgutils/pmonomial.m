function [f,fx,fy,fz]=pmonomial(x,p)
%PMONOMIAL  Multivariate monomials and derivatives
%
%    Compute all multivariate monomials prod(x_i^p_i) up to the
%    maximum degree sum(p_i) = p for given x_i.
%
%    Syntax: [f,fx,fy,fz]=pmonomial(x,p)
%
%    x:  coordinates for evaluation, dimensions [npoints,dim]
%    p:  maximum monomial order
%    f:  values, dimensions [npoints,nmonomials]
%    fs: s-derivatives, dimensions [npoints,nmonomials]

dim=size(x,2);

pp=prod((p+(1:dim))./(1:dim));
f=zeros(size(x,1),pp);
fx=zeros(size(x,1),pp);
fy=zeros(size(x,1),pp);
fz=zeros(size(x,1),pp);
pos=1;

switch dim
 case 0
  f=0*f+1;
 case 1
  for i=0:p
    f(:,pos)=x(:,1).^i;
    if i>=1
      fx(:,pos)=i*x(:,1).^(i-1);
    end
    pos=pos+1;
  end
 case 2
  for i=0:p
    for j=0:i
      f(:,pos)=x(:,1).^(i-j).*x(:,2).^j;
      if i-j>=1
        fx(:,pos)=(i-j)*x(:,1).^(i-j-1).*x(:,2).^j;
      end
      if j>=1
        fy(:,pos)=j*x(:,1).^(i-j).*x(:,2).^(j-1);
      end
      pos=pos+1;
    end
  end
 case 3
  for i=0:p
    for j=0:p
      for k=0:p
        if i+j+k>=0 & i+j+k<=p
          f(:,pos)=x(:,1).^i.*x(:,2).^j.*x(:,3).^k;
          if i>=1
            fx(:,pos)=i*x(:,1).^(i-1).*x(:,2).^j.*x(:,3).^k;
          end
          if j>=1
            fy(:,pos)=j*x(:,1).^i.*x(:,2).^(j-1).*x(:,3).^k;
          end
          if k>=1
            fz(:,pos)=k*x(:,1).^i.*x(:,2).^j.*x(:,3).^(k-1);
          end
          pos=pos+1;
        end
      end
    end
  end
 otherwise
   error('Dimension not implemented.');
end
