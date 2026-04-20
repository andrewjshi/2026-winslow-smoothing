function [s,tlocal,sbnd,tbndlocal]=lagrangepnts(p,dim,eltype,s0)

if nargin<2, dim=2; end
if nargin<3, eltype=t_simplex; end

if eltype==t_simplex
    if nargin >= 4 
        if norm(s0 - (0:p)'/p, inf) > 1e-12
            error('Only equidistant nodes supported for simplex elements');
        end
    end
    if p==0
        s=ones(1,dim+1)/(dim+1);
        tlocal=ones(1,dim+1);
    else
        p=double(p);
        switch dim
          case 0
            s=1;
            tlocal=1;
          case 1
            s=[(0:p)',(p:-1:0)']/p;
            if nargout>=2
                tlocal=[1:p;2:p+1]';
            end
          case 2
            [u,v]=ndgrid((0:p)/p,(0:p)/p);
            s=[u(:),v(:)];
            s=[s,1-sum(s,2)];
            s=s(s(:,3)>=0,:);
            if nargout>=2
                tlocal=zeros(0,3);
                loc=0;
                for ii=1:p
                    jj=p+1-ii;
                    t1=[1:jj; 2:jj+1; jj+1+(1:jj)]';
                    t2=[2:jj; jj+1+(2:jj); jj+1+(1:jj-1)]';
                    tlocal=[tlocal; t1+loc; t2+loc];
                    loc=loc+jj+1;
                end
            end
          case 3
            [u,v,w]=ndgrid((0:p)/p,(0:p)/p,(0:p)/p);
            s=[u(:),v(:),w(:)];
            s=[s,1-sum(s,2)];
            s=s(s(:,4)>=0,:);
            if nargout>=2
                [foo,foo,tlocal]=tetrahedronmesh(p+1);
            end
          otherwise
            error('Dimension not implemented.');
        end
    end
    if nargout>=3
        [sbnd,tbndlocal]=lagrangepnts(p,dim-1,eltype);
    end
elseif eltype==t_block
    if nargin < 4
        s0 = (lobattopoints(p+1)+1)/2;
    end
    if p==0
        s=ones(1,2^dim)/2^dim;
        tlocal=ones(1,2^dim);
    else
        switch dim
          case 0
            s=1;
            tlocal=1;
          case 1
            s=s0(:);
            s=[1-s, s];
            if nargout>=2
                tlocal=[1:p;2:p+1]';
            end
          case 2
            [s,t]=ndgrid(s0,s0);
            s=s(:); t=t(:);
            s=[(1-s).*(1-t), ...
               (  s).*(1-t), ...
               (1-s).*(  t), ...
               (  s).*(  t)];
            if nargout>=2
                ix=reshape(1:(p+1)^2,(p+1)*ones(1,2));
                tlocal=[ix(1:p,1:p),ix(2:p+1,1:p),ix(1:p,2:p+1),ix(2:p+1,2:p+1)];
                tlocal=reshape(tlocal,[],4);
            end
          case 3
            [s,t,u]=ndgrid(s0,s0,s0);
            s=s(:); t=t(:); u=u(:);
            s=[(1-s).*(1-t).*(1-u), ...
               (  s).*(1-t).*(1-u), ...
               (1-s).*(  t).*(1-u), ...
               (  s).*(  t).*(1-u), ...
               (1-s).*(1-t).*(  u), ...
               (  s).*(1-t).*(  u), ...
               (1-s).*(  t).*(  u), ...
               (  s).*(  t).*(  u)];
            if nargout>=2
                ix=reshape(1:(p+1)^3,(p+1)*ones(1,3));
                tlocal=cat(3, ...
                           ix(1:p,1:p,1:p), ...
                           ix(2:p+1,1:p,1:p), ...
                           ix(1:p,2:p+1,1:p), ...
                           ix(2:p+1,2:p+1,1:p), ...
                           ix(1:p,1:p,2:p+1), ...
                           ix(2:p+1,1:p,2:p+1), ...
                           ix(1:p,2:p+1,2:p+1), ...
                           ix(2:p+1,2:p+1,2:p+1));
                tlocal=reshape(tlocal,[],8);
            end
          otherwise
            error('Dimension not implemented.');
        end
    end
    if nargout>=3
        if nargin >= 4
            [sbnd,tbndlocal]=lagrangepnts(p,dim-1,eltype,s0);
        else
            [sbnd,tbndlocal]=lagrangepnts(p,dim-1,eltype);
        end
    end
else
    error('Unknown element type');
end
 
if nargout>=2
    tlocal=int32(tlocal);
end
