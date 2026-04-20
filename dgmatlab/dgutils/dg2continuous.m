function u = dg2continuous(msh, u, method, solver)

if nargin<3, method = 'average'; end
if nargin<4, solver = '\'; end

switch lower(method)
  case 'average'
    if isstruct(msh)
        p1 = msh.p1;
    else
        p1 = msh;
    end
    
    [q,W] = mapdg2cg(p1);
    [iW,jW,sW] = find(W);
    W0 = sparse(iW, jW, 1, size(W,1), size(W,2));
    
    u = permute(reshape(W0' * (W * reshape(permute(u, [1,3,2]), size(W,2), ...
                            [])), size(u,1), size(u,3), size(u,2)), [1,3,2]);
  case 'l2'
    data = dginit(msh);
    [ns,nc,nt] = size(u);
    phi = squeeze(data.gfs(:,1,:))';
    wphi = squeeze(data.wgfs(:,1,:))';

    J = geojac(msh, data);
    ug = phi * reshape(permute(u, [1,3,2]), [ns,nt*nc]);

    res = wphi' * (repmat(J, [1,nc]) .* ug) / factorial(size(msh.p,1));
    [ix, ~, W1] = mapdg2cg(msh.p1);
    ix = repmat(ix(:), 1, nc);
    comp = repmat(1:nc, ns*nt, 1);

    MD = dgmass(msh, data);
    [Dii,Djj] = dgindices(msh, data, 'mass', 1);
    CM = sparse(ix(Dii), ix(Djj), MD);
    rhs = full(sparse(ix(:), comp(:), res(:), max(ix(:)), nc));
    switch solver
      case '\'
        uc = CM \ rhs;
      case 'cg'
        %RCM = ichol(CM);
        RCM = ichol(CM, struct('type', 'ict', 'droptol', 1e-3));
        uc = 0 * rhs;
        for ic = 1:nc
            [uc(:,ic), flag] = pcg(CM, rhs(:,ic), 1e-8, 1000, RCM, RCM');
            if flag ~= 0
                error('No convergence in CG');
            end
        end
      otherwise
        error('Unknown solver.');
    end

    u = permute(reshape(W1'*uc, ns, nt, nc), [1,3,2]);
  otherwise
    error('Unknown method');
end
