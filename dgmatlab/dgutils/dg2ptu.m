function [p0,t0,u0] = dg2ptu(msh, u, utype, nref)
%DG2PTU  High-order 3DG msh + solution to linear p,t,u format
%
%   [p0,t0,u0] = dg2ptu(msh, u, utype, nref)

    if nargin < 3, utype = 1; end
    if nargin < 4, nref = 0; end

    [nv,nt] = size(msh.t);
    dim = size(msh.p,1);
    ns = size(msh.s,1);
    porder = double(msh.porder);

    if porder == 0
        [ss,tt] = lagrangepnts(1,dim,msh.eltype);
        pp = zeros(nv,dim,nt);
        for ic = 1:dim
            pp(:,ic,:) = reshape(msh.p(ic,msh.t + 1), nv, 1, nt);
        end
        u = repmat(u, [nv,1,1]);
    else
        ss = msh.s;
        tt = msh.tlocal;
        if min(tt(:)) == 0
            tt = tt + 1;
        end
        pp = msh.p1;
    end
    
    if size(u,1) == 1
        u = repmat(u, [ns,1,1]);
    end

    u = dgeval(u, utype, dim);

    switch msh.eltype
      case t_simplex
        if nref > 0
            A0 = pmonomial(ss(:,1:dim), porder);
            [ss,tt] = lagrangepnts(porder * 2^nref, dim);
            A = pmonomial(ss(:,1:dim), porder) / A0;
            
            nss = size(ss,1);
            sz = size(pp); if length(sz) == 2, sz = [sz,1]; end
            pp = reshape(A * reshape(pp, ns, sz(2)*sz(3)), [nss,sz(2),sz(3)]);
            sz = size(u); if length(sz) == 2, sz = [sz,1]; end
            u = reshape(A * reshape(u, ns, sz(2)*sz(3)), [nss,sz(2),sz(3)]);
        end
    
        nss = size(ss,1);
        p0 = reshape(permute(pp,[1,3,2]), [ nss*nt,dim]);
        t0 = kron(ones(nt,1,'int32'), tt) + kron(nss*int32(0:nt-1)', 0*tt + 1);
        u0 = squeeze(u);

      case t_block
        ninterp = porder * 2^nref + 1;
        V = plegendre(linspace(-1,1,ninterp)',porder) / plegendre(2*msh.s0 - 1, porder);
        error('TODO: General dim');
        for it = 1:nt
            x = V*reshape(msh.p1(:,1,it),[porder+1,porder+1])*V';
            y = V*reshape(msh.p1(:,2,it),[porder+1,porder+1])*V';
            c = V*reshape(u(:,:,it),[porder+1,porder+1])*V';
        end
      otherwise
        error('Unknown element type');
    end
end

