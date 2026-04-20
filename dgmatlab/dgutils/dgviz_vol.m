function viz = dgviz_vol(msh, u, navstok)
%DGVIZ_VOL Create volumetric vizualization data
%
%   viz = dgviz_vol(msh, u=[], navstok=0)
%
%   See DGVTKWRITE for examples.

    if nargin < 2, u = []; end
    if nargin < 3, navstok = 0; end

    np=size(msh.p,2);
    nt=size(msh.t,2);
    ns=size(msh.s,1);
    dim=size(msh.p,1);
    
    porder=double(msh.porder);
    
    if msh.eltype == t_block
        % Interpolate to equidistant points
        s0 = (0:porder) / porder;
        msh.p1 = dginterp(msh.p1, msh.s0, s0, dim, msh.eltype);
        if ~isempty(u)
            u = dginterp(u, msh.s0, s0, dim, msh.eltype);
        end
        msh.s0 = s0;
    end
    
    p1 = reshape(permute(msh.p1,[2,1,3]),dim,[])';
    t1 = reshape(1:ns*nt,ns,nt)';
    
    if ~isempty(u)
        u = reshape(permute(u,[2,1,3]),size(u,2),[])';
    
        if navstok >= 1
            % Navier-Stokes
            u1 = {{'scalars', 'Density_r', u(:,1)}, ...
                  {'vectors', 'Momentum_ru', u(:,2:dim+1)}};
            if navstok >= 2 & size(u,2) >= dim+2
                u1{3} = {'scalars', 'Energy_rE', u(:,dim+2)};
            end
        else
            u1 = u;
        end     
    end
    
    viz = struct('porder',porder, 'eltype',msh.eltype, 'dim', dim, ...
                 'p',p1, 't',t1);
    if ~isempty(u)
        viz.u = u1;
    end
end
