% w2_3d.m - 3D analogue of w2.m.
%
% Structured tetrahedral reference mesh of the unit cube,
% violent perturbation on face 4 (y=1), Winslow-Picard smoothing to
% untangle. Saves three figures:
%   figures/w2_3d/mesh_reference.png
%   figures/w2_3d/mesh_tangled.png
%   figures/w2_3d/mesh_smoothed.png

porder = 4;
n = 6;
msh = mshcube(n+1, n+1, n+1);
msh = nodealloc(msh, porder);

figdir = 'figures/w2_3d';
if ~exist(figdir, 'dir'), mkdir(figdir); end

ix_all = find(any(ismember(msh.t2t, -(1:6)), 1));

% --- Figure 1: reference mesh -------------------------------------------
figure(1); clf;
tetplotcurved(msh, 1:6, ix_all, [.8, 1, .8]);
delete(findall(gcf, 'type', 'colorbar'));
plotbndnodes(msh, 1:6);
axis equal; axis off; view(-130, 20); drawnow;
set(gcf, 'Position', [114 1 560 420]);
exportgraphics(gcf, fullfile(figdir, 'mesh_reference.png'), 'Resolution', 200);

% --- Build the perturbation: displace nodes on boundary 1 (x=0) ----------
% Matches the 2D w2 pattern: push the face inward into the domain.
% Grid spacing is Delta = 1/n; amplitude cuts through several element
% layers producing heavy tangling near the perturbed face.
[~, edg] = getbndnodes(msh, dginit(msh), 1);
x = msh.p1(:,1,:);
y = msh.p1(:,2,:);
z = msh.p1(:,3,:);
newx = x;
newy = y;
newz = z;
newx(edg) = x(edg) + 0.4*sin(pi*y(edg)).*sin(pi*z(edg));
curvep1 = cat(2, newx, newy, newz);

msh_tangled = msh;
msh_tangled.p1 = curvep1;

% --- Figure 2: tangled initial configuration -----------------------------
figure(2); clf;
tetplotcurved(msh_tangled, 1:6, ix_all, [.8, 1, .8]);
delete(findall(gcf, 'type', 'colorbar'));
plotbndnodes(msh_tangled, 1:6);
axis equal; axis off; view(-130, 20); drawnow;
set(gcf, 'Position', [114 1 560 420]);
exportgraphics(gcf, fullfile(figdir, 'mesh_tangled.png'), 'Resolution', 200);

% --- Run Winslow-Picard ---------------------------------------------------
doplot = false;
msh1 = elliptic_smoothing(msh, curvep1, doplot);

% --- Figure 3: final smoothed mesh ---------------------------------------
figure(3); clf;
tetplotcurved(msh1, 1:6, ix_all, [.8, 1, .8]);
delete(findall(gcf, 'type', 'colorbar'));
plotbndnodes(msh1, 1:6);
axis equal; axis off; view(-130, 20); drawnow;
set(gcf, 'Position', [114 1 560 420]);
exportgraphics(gcf, fullfile(figdir, 'mesh_smoothed.png'), 'Resolution', 200);

% --- Scaled Jacobian histograms (Fortunato-Persson 2016, Section 4) -----
% Per element I_e = min J(xi) / max J(xi), approximated at Gauss points.
% I_e <= 1; straight-sided simplex => I_e = 1; I_e <= 0 => inverted.
I_ref      = min(geojac(msh))       ./ max(geojac(msh));
I_tangled  = min(geojac(msh_tangled)) ./ max(geojac(msh_tangled));
I_smoothed = min(geojac(msh1))      ./ max(geojac(msh1));

Iall = [I_ref(:); I_tangled(:); I_smoothed(:)];
edges = linspace(min(Iall), 1, 40);

configs = {
    'reference', I_ref,      'hist_reference.png';
    'tangled',   I_tangled,  'hist_tangled.png';
    'smoothed',  I_smoothed, 'hist_smoothed.png';
};
for k = 1:3
    figure(3+k); clf;
    histogram(configs{k,2}(:), edges, 'FaceColor', [.6, .85, .6]);
    xline(0, 'k--');
    xlabel('scaled Jacobian'); ylabel('element count');
    title(sprintf('%s: min I = %.4g', configs{k,1}, min(configs{k,2})));
    set(gcf, 'Position', [114 1 560 420]); drawnow;
    exportgraphics(gcf, fullfile(figdir, configs{k,3}), 'Resolution', 200);
end

fprintf('\nSaved figures to %s/\n', figdir);

function plotbndnodes(msh, bnds)
data = dginit(msh);
xyz = [];
for it = 1:size(msh.t, 2)
    for j = 1:4
        if ismember(-msh.t2t(j, it), bnds)
            cegix = data.egix(:, j) + 1;
            xyz = [xyz; squeeze(msh.p1(cegix, :, it))]; %#ok<AGROW>
        end
    end
end
xyz = unique(xyz, 'rows');
hold on;
plot3(xyz(:,1), xyz(:,2), xyz(:,3), '.', 'Color', [0, 0.3, 0.9], 'MarkerSize', 5);
hold off;
end
