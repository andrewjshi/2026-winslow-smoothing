function [phi_map, phi_ix] = linelast_deform(msh, bndfix, bndx, bndy, fphi)
%LINELAST_DEFORM
%    [PHI_MAP, PHI_IX] = LINELAST_DEFORM(MSH, BNDFIX, BNDX, BNDY, FPHI)
%    [PHI_MAP, PHI_IX] = LINELAST_DEFORM(MSH, IXFIX, IXMU)
%
%    Example (cylmsh1):
%       msh = h5freadstruct('cylmsh1p3.h5');
%
%       fphi = @(p) abs(dsegment(p, [-1,5; -2.5,0])) < 1e-3;
%       if 0 % Old syntax
%         bndfix = 2; bndx = [1,4]; bndy = [3,5];
%         [phi_map, phi_ix] = linelast_deform(msh, bndfix, bndx, bndy, fphi);
%       else % New syntax
%         xy = dg2cg(msh, msh.p1);
%         ncg = size(xy,1);
%         e = getbndnodes(msh, dginit(msh), 2);
%         ex = getbndnodes(msh, dginit(msh), [1,4]);
%         ey = getbndnodes(msh, dginit(msh), [3,5]);
%         ixfix = unique([e; e+ncg; ex+ncg; ey]);
%         emu = find(fphi(xy));
%         ixmu = setdiff([emu; emu+ncg], ixfix);
%         [phi_map, phi_ix] = linelast_deform(msh, ixfix, ixmu);
%       end
%       
%       phi = 0.05 * randn(size(phi_map,2),1);  % random
%       phi(1:end/2+1) = phi(1:end/2+1) + 0.5;  % x shift
%       dx = phi_map * phi;
%
%       p1cg = dg2cg(msh, msh.p1);
%       p1cg(phi_ix) = p1cg(phi_ix) + dx;
%       newp1 = cg2dg(msh, p1cg);
%
%       msh1 = msh; msh1.p1 = newp1; dgmeshplot_curved(msh1, 2, 1)

if nargin == 3
    ixfix = bndfix;
    ixmu = bndx;
    syntax = 1; % New syntax
elseif nargin == 5
    syntax = 0; % Old syntax
else
    error('Incorrect syntax');
end

dim = size(msh.p,1);
porder = double(msh.porder);
data = dginit(msh);

xycg = dg2cg(msh, msh.p1);
ncg = size(xycg,1);

% Find fixed and mu DOFs

switch syntax
  case 0  % Old syntax
    e = getbndnodes(msh, data, bndfix);
    ex = getbndnodes(msh, data, bndx);
    ey = getbndnodes(msh, data, bndy);
    
    e1 = e;
    for idim = 1:dim-1
        e1 = [e1; e + idim*ncg];
    end
    e1 = [e1; ex + ncg];
    e1 = [e1; ey];
  
    if isnumeric(fphi)
        phiix = getbndnodes(msh, data, fphi);
    else
        phiix = find(fphi(xycg));
    end
    phiix = [phiix; phiix + ncg];
    phiix = setdiff(phiix, e1);
  case 1  % New syntax
    
    e1 = ixfix;
    phiix = ixmu;
    
    if ~isempty(intersect(e1,phiix))
        error('Overlapping indices in ixfix and ixmu');
    end
end

ixkeep = setdiff((1:dim*ncg)', e1);
[~,phiixkeep] = ismember(phiix, ixkeep);

% Create elimination matrices
[r,W,W1] = mapdg2cg(msh.p1);
Kperm = mkKperm(msh);
[foo,Kperminv] = sort(Kperm);
W11 = kron(speye(dim,dim), W1);
W11perm = W11(:,Kperminv);
W11permkeep = W11perm(ixkeep,:);
[Dii,Djj] = dgindices(msh, data, 'blockdiag', dim);

% Linear elasticity physics parameters
E      = 1e0;
nu     = 0;
lambda = E*nu/(1+nu)/(1-2*nu);
mu     = E/2/(1+nu);
rho    = 1.0;
clear pars
phys.pars = [lambda,mu,rho];
phys.bndcnds = int32(ones(1,100));

% Assemble and eliminate
[R,DA] = dglinelast(0*msh.p1, msh, data, phys, 'cg');
F = W11permkeep * R(:);
K = sparse(Dii(:), Djj(:), DA(:));
K = W11permkeep * K * W11permkeep';

% Form inverse mapping
nphiix = numel(phiix);
e = 1:nphiix;
nkeep = numel(ixkeep);
B = sparse(e, phiixkeep, 1, nphiix, nkeep);

phi_map = full(K \ (B' / (B*(K\B'))));
phi_ix = ixkeep;
