% w_sliver_only.m - W-modification applied ONLY to sliver-cluster elements.
% Everywhere else W_K = I, so those elements behave like unmodified Winslow.
% Tests whether the sliver correction is being drowned out by averaging.

porder = 4;
n = 10;
msh = mshsquare(n+1, n+1);
msh = nodealloc(msh, porder);
moves = [0.8 0.8 0.72 0.79];
[~, ix] = min(vecnorm(double(msh.p) - moves(1:2)', 2, 1));
msh.p(:, ix) = moves(3:4);
msh = nodealloc(msh, porder);

[~, edg] = getbndnodes(msh, dginit(msh), 4);
xx = msh.p1(:,1,:); yy = msh.p1(:,2,:);
newx = xx; newy = yy;
newx(edg) = xx(edg) + 0.5*sin(pi*yy(edg));
curvep1 = cat(2, newx, newy);

% Identify sliver elements by eta > threshold in the reference
data = dginit(msh);
p1x = permute(msh.p1(:,1,:), [1,3,2]);
p1y = permute(msh.p1(:,2,:), [1,3,2]);
phiX = permute(data.gfs(:,2,:), [1,3,2]);
phiY = permute(data.gfs(:,3,:), [1,3,2]);
xX = phiX'*p1x; xY = phiY'*p1x;
yX = phiX'*p1y; yY = phiY'*p1y;
detJ = xX.*yY - xY.*yX;
Jfro2 = xX.^2 + xY.^2 + yX.^2 + yY.^2;
etap = Jfro2 ./ (2 * max(abs(detJ), eps));
eta_elem = max(etap, [], 1)';
sliver_mask = eta_elem > 2.0;
fprintf('sliver cluster: %d / %d elements (eta > 2)\n', sum(sliver_mask), numel(sliver_mask));

% Build W_K: identity everywhere, equilateral target on sliver elements only.
d = 2; nt = size(msh.t, 2); ns = size(msh.s, 1);
WK_full = compute_WK_simplex(msh, 'equilateral');
WK_sparse = repmat(eye(d), [1, 1, nt]);
WK_sparse(:, :, sliver_mask) = WK_full(:, :, sliver_mask);

figdir = 'figures/sliver_only'; if ~exist(figdir,'dir'), mkdir(figdir); end

for mode = {'const', 'p1', 'vhp'}
    m = mode{1};
    switch m
    case 'const'
        W_h = zeros(ns, d, d, nt);
        for it = 1:nt
            for i = 1:ns
                W_h(i,:,:,it) = WK_sparse(:,:,it);
            end
        end
    case {'p1','vhp'}
        msh_tmp = msh;
        msh_tmp.W_K_sparse = WK_sparse;
        W_h = project_sparse(msh, WK_sparse, m);
    end
    msh.W_h = W_h;
    try
        msh1 = elliptic_smoothing_W(msh, curvep1, false);
        [~,ix_sliver] = min(vecnorm(double(msh.p) - [0.72; 0.79], 2, 1));
        pflat = reshape(permute(msh1.p1,[1,3,2]), [], 2);
        pref  = reshape(permute(msh.p1,[1,3,2]), [], 2);
        maskS = vecnorm(pref - [0.72 0.79], 2, 2) < 1e-8;
        idxS = find(maskS);
        xS = mean(pflat(idxS,1)); yS = mean(pflat(idxS,2));

        data1 = dginit(msh1);
        p1x1 = permute(msh1.p1(:,1,:),[1,3,2]); p1y1 = permute(msh1.p1(:,2,:),[1,3,2]);
        xX1=phiX'*p1x1; xY1=phiY'*p1x1; yX1=phiX'*p1y1; yY1=phiY'*p1y1;
        detJ1 = xX1.*yY1-xY1.*yX1;
        Jfro21 = xX1.^2+xY1.^2+yX1.^2+yY1.^2;
        etap1 = Jfro21 ./ (2*max(abs(detJ1),eps));
        eta1 = max(etap1,[],1)';
        fprintf('[%-5s] sliver image = (%.4f, %.4f), max eta = %.3g, mean eta = %.3g\n', ...
            m, xS, yS, max(eta1), mean(eta1));

        figure; clf; dgmeshplot_curved(msh1, 4, 0, 0);
        title(sprintf('sliver-only W (%s)', m));
        set(gcf,'Position',[100 100 600 600]);
        exportgraphics(gcf, fullfile(figdir, sprintf('mesh_%s.png', m)), 'Resolution',200);
    catch ME
        fprintf('[%-5s] FAILED: %s\n', m, ME.message);
    end
end

function W_h = project_sparse(msh, WK, kind)
d = size(msh.p,1); nt = size(msh.t,2); ns = size(msh.s,1);
t = double(msh.t)+1;
absK = zeros(1,nt);
for it=1:nt
    corners = msh.p(:, t(:,it));
    E = corners(:, 2:end) - corners(:, 1);
    absK(it) = abs(det(E))/factorial(d);
end
if strcmp(kind,'p1')
    npv = size(msh.p,2);
    numW = zeros(d,d,npv); denW = zeros(1,npv);
    for it=1:nt
        for k=1:size(t,1)
            vi=t(k,it);
            numW(:,:,vi)=numW(:,:,vi)+absK(it)*WK(:,:,it);
            denW(vi)=denW(vi)+absK(it);
        end
    end
    W_CG = numW ./ reshape(denW,1,1,[]);
    S = msh.s(:,1:d); phi_lin = [1-sum(S,2), S];
    W_h = zeros(ns,d,d,nt);
    for it=1:nt
        ci = t(:,it);
        Wc = W_CG(:,:,ci);
        for i=1:ns
            acc = zeros(d,d);
            for k=1:(d+1), acc = acc + phi_lin(i,k)*Wc(:,:,k); end
            W_h(i,:,:,it) = acc;
        end
    end
else  % vhp
    q = mapdg2cg(msh.p1);
    ncg = max(q);
    numW = zeros(d,d,ncg); denW = zeros(1,ncg);
    for it=1:nt
        w = absK(it)*WK(:,:,it);
        for i=1:ns
            cg = q(i+ns*(it-1));
            numW(:,:,cg) = numW(:,:,cg) + w;
            denW(cg) = denW(cg) + absK(it);
        end
    end
    W_CG = numW ./ reshape(denW,1,1,[]);
    W_h = zeros(ns,d,d,nt);
    for it=1:nt
        for i=1:ns
            cg = q(i+ns*(it-1));
            W_h(i,:,:,it) = W_CG(:,:,cg);
        end
    end
end
end
