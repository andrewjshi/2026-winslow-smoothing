function viz = dgviz_bnd(msh, u, app, phys, postqtys, bndnbrs)
%DGVIZ_BND Create boundary vizualization data (needs application data)
%
%   viz = dgviz_bnd(msh, u, app, phys, postqtys, bndnbrs='all')
%
%   See DGVTKWRITE for examples.

    if nargin < 6, bndnbrs = 'all'; end

    np=size(msh.p,2);
    nt=size(msh.t,2);
    nes=size(msh.sbnd,1);
    dim=size(msh.p,1);
    
    porder=double(msh.porder);
    
    if strcmp(bndnbrs, 'all')
        bndnbrs = unique(-msh.t2t(msh.t2t<0));
    end
        
    data = dginit(msh, struct('gaussadds',true));
    p1 = [];
    u1 = [];
    for ibnd = bndnbrs(:)'
        xy = app(u, msh, data, phys, ['posteval x ', int2str(ibnd)]);
        p1 = cat(3, p1, xy);
        cu = [];
        for postqty = postqtys(:)'
            cu1 = app(u, msh, data, phys, ['posteval ', int2str(postqty), ' ', int2str(ibnd)]);
            cu = cat(2, cu, cu1);
        end
        u1 = cat(3, u1, cu);
    end
    p1 = p1(end-nes+1:end,:,:);
    u1 = u1(end-nes+1:end,:,:);
    ne = size(p1,3);
    
    if msh.eltype == t_block
        % TODO: Need to implement gaussadds in dginit for hexes (tensor-based)
        error('Not implemented for hexes');
        
        % Interpolate to equidistant points
        s0 = (0:porder) / porder;
        p1 = dginterp(p1, msh.s0, s0, dim-1, msh.eltype);
        u1 = dginterp(u1, msh.s0, s0, dim-1, msh.eltype);
        msh.s0 = s0;
    end
    
    p1 = reshape(permute(p1,[2,1,3]),dim,[])';
    t1 = reshape(1:nes*ne,nes,ne)';
    
    u1 = reshape(permute(u1,[2,1,3]),size(u1,2),[])';
    
    viz = struct('porder',porder, 'eltype',msh.eltype, 'dim',dim-1, ...
                 'p',p1, 't',t1);
    viz.u = u1;
end
