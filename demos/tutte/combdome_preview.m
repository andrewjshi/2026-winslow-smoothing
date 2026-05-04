% combdome_preview.m -- Preview the dome reference shape and the
% comb -> dome boundary correspondence. No Tutte, no smoothing, no
% Winslow. Just pictures.
%
% Dome = comb base rectangle [0, W] x [0, H_base] capped by a half-
% ellipse on top. Major semi-axis a = W/2 (locked by rectangle width).
% Minor semi-axis b chosen so the half-ellipse arc length equals the
% comb's teeth-region arc length. Mapping is isometric by arc length.

% --- Comb parameters (must match comb.m) --------------------------------
W       = 0.92;
H_base  = 0.20;
H_top   = 0.38;
margin  = 0.10;
tooth_w = 0.18;
gap_w   = 0.10;
n_teeth = 3;

% --- Comb corners (CCW) -------------------------------------------------
C = [0, 0; W, 0; W, H_base];
x = W - margin;
for k = 1:n_teeth
    C = [C; x, H_base; x, H_top; x - tooth_w, H_top; x - tooth_w, H_base]; %#ok<AGROW>
    x = x - tooth_w - gap_w;
end
C = [C; 0, H_base];
nc = size(C, 1);

% --- Edge lengths and teeth-region arc length ---------------------------
edge_vec = C([2:nc 1], :) - C;
edge_len = vecnorm(edge_vec, 2, 2);

% Edge i goes from corner i to corner mod(i,nc)+1.
% Bottom = edge 1, right = edge 2, teeth-region = edges 3..nc-1, left = edge nc.
L_bottom = edge_len(1);
L_right  = edge_len(2);
L_left   = edge_len(nc);
L_teeth  = sum(edge_len(3:nc-1));
fprintf('comb: bottom=%.3f, right=%.3f, left=%.3f, teeth=%.3f\n', ...
    L_bottom, L_right, L_left, L_teeth);

% --- Dome geometry ------------------------------------------------------
a_dome = W/2;
b_dome = solve_b_for_arc(L_teeth, a_dome);
fprintf('dome: a=%.4f, b=%.4f, half-perim=%.4f (target %.4f)\n', ...
    a_dome, b_dome, ellipse_half_perim(a_dome, b_dome), L_teeth);

% --- Dense theta grid for cap, with cumulative arc length ---------------
N_theta = 4000;
theta_grid = linspace(0, pi, N_theta+1)';
ds_dtheta  = sqrt(a_dome^2 * sin(theta_grid).^2 + b_dome^2 * cos(theta_grid).^2);
s_grid     = [0; cumsum((ds_dtheta(1:end-1) + ds_dtheta(2:end))/2 * (pi/N_theta))];
% By construction, s_grid(end) ~= L_teeth (Ramanujan vs. trapezoidal).
% Renormalize so the endpoints coincide exactly.
s_grid = s_grid * (L_teeth / s_grid(end));

% --- Map every comb corner onto the dome (isometric by arc length) ------
C_dome = zeros(nc, 2);
for i = 1:nc
    if i == 1 || i == 2 || i == 3
        % Corners 1, 2, 3 = (0,0), (W,0), (W,H_base) -- match dome
        C_dome(i, :) = C(i, :);
    elseif i == nc
        % Corner nc = (0, H_base) -- matches dome's top-left of rectangle
        C_dome(i, :) = C(i, :);
    else
        % Teeth-region corners: arc length from corner 3
        s = sum(edge_len(3:i-1));
        theta = interp1(s_grid, theta_grid, s, 'linear');
        C_dome(i, :) = [W/2 + a_dome*cos(theta), H_base + b_dome*sin(theta)];
    end
end

% --- Sample boundary points along each edge for a smooth picture --------
n_per_edge_density = 25;
bnd_comb = [];
bnd_dome = [];
for i = 1:nc
    L = edge_len(i);
    npts = max(2, round(L * n_per_edge_density));
    for k = 0:npts-1
        t = k / npts;
        bnd_comb = [bnd_comb; C(i,:) + t*edge_vec(i,:)]; %#ok<AGROW>
        if i == 1 || i == 2 || i == nc
            bnd_dome = [bnd_dome; C(i,:) + t*edge_vec(i,:)]; %#ok<AGROW>
        else
            % Teeth-region edge: arc length from corner 3 to this point
            s = sum(edge_len(3:i-1)) + t*L;
            theta = interp1(s_grid, theta_grid, s, 'linear');
            bnd_dome = [bnd_dome; [W/2+a_dome*cos(theta), H_base+b_dome*sin(theta)]]; %#ok<AGROW>
        end
    end
end
bnd_comb = [bnd_comb; bnd_comb(1,:)];   % close
bnd_dome = [bnd_dome; bnd_dome(1,:)];

% --- Color-code each tooth and gap --------------------------------------
% Edge groups for the teeth region:
%   margin        = edge 3
%   tooth k       = edges 3 + 1 + (k-1)*4 .. 3 + 3 + (k-1)*4   (3 edges per tooth)
%   gap k         = edge 3 + 4 + (k-1)*4   (between teeth)
%   final stub    = edge nc-1
% Easier: classify each edge as 'rect-bottom', 'rect-right', 'rect-left',
% 'margin', 'tooth k', 'gap k', 'stub'.
edge_class = strings(nc, 1);
edge_class(1) = "rect_bottom";
edge_class(2) = "rect_right";
edge_class(nc) = "rect_left";
edge_class(3) = "margin";
for k = 1:n_teeth
    base = 3 + (k-1)*4;
    edge_class(base+1) = sprintf("tooth%d", k);
    edge_class(base+2) = sprintf("tooth%d", k);
    edge_class(base+3) = sprintf("tooth%d", k);
    if k < n_teeth
        edge_class(base+4) = sprintf("gap%d", k);
    end
end
edge_class(nc-1) = "stub";

% --- Plot ---------------------------------------------------------------
figdir = 'figures/combdome';
if ~exist(figdir, 'dir'), mkdir(figdir); end

f = figure('Position', [50 50 1300 600]); clf;

tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% Tooth colors
tooth_colors = lines(n_teeth);

% Left: comb with color-coded edges
ax1 = nexttile; hold on; axis equal; axis off;
title('comb boundary');
plot_edges(C, edge_class, tooth_colors, n_teeth);
plot(C(:,1), C(:,2), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
% Label corners
for i = 1:nc
    text(C(i,1)+0.005, C(i,2)+0.012, sprintf('%d', i), 'FontSize', 8);
end

% Right: dome with same color-coded edges
ax2 = nexttile; hold on; axis equal; axis off;
title(sprintf('dome reference (a=%.3f, b=%.3f)', a_dome, b_dome));
% Plot dome boundary as fine sampled polylines, color by edge class
% We sample each comb edge with many points and connect through C_dome
for i = 1:nc
    npts = 60;
    if i == 1 || i == 2 || i == nc
        % Straight segment
        p0 = C_dome(i,:); p1 = C_dome(mod(i,nc)+1,:);
        ts = linspace(0, 1, npts)';
        seg = (1-ts).*p0 + ts.*p1;
    else
        % Arc: walk arc length from corner i to corner i+1
        s0 = sum(edge_len(3:i-1));
        s1 = s0 + edge_len(i);
        ss = linspace(s0, s1, npts)';
        ths = interp1(s_grid, theta_grid, ss, 'linear');
        seg = [W/2 + a_dome*cos(ths), H_base + b_dome*sin(ths)];
    end
    plot(seg(:,1), seg(:,2), '-', 'Color', edge_color(edge_class(i), tooth_colors, n_teeth), 'LineWidth', 2);
end
plot(C_dome(:,1), C_dome(:,2), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
for i = 1:nc
    text(C_dome(i,1)+0.012, C_dome(i,2)+0.012, sprintf('%d', i), 'FontSize', 8);
end

linkaxes([ax1 ax2], 'off');
axis(ax1, [-0.05, W+0.05, -0.05, max(H_top, H_base+b_dome)+0.05]);
axis(ax2, [-0.05, W+0.05, -0.05, max(H_top, H_base+b_dome)+0.05]);

exportgraphics(gcf, fullfile(figdir, 'combdome_preview.png'), 'Resolution', 200);
fprintf('\nSaved %s\n', fullfile(figdir, 'combdome_preview.png'));


% =========================================================================
function plot_edges(C, edge_class, tooth_colors, n_teeth)
nc = size(C, 1);
for i = 1:nc
    p0 = C(i, :); p1 = C(mod(i,nc)+1, :);
    plot([p0(1) p1(1)], [p0(2) p1(2)], '-', ...
        'Color', edge_color(edge_class(i), tooth_colors, n_teeth), ...
        'LineWidth', 2);
end
end

function c = edge_color(cls, tooth_colors, n_teeth) %#ok<INUSD>
switch char(cls)
    case 'rect_bottom', c = [0.4 0.4 0.4];
    case 'rect_right',  c = [0.4 0.4 0.4];
    case 'rect_left',   c = [0.4 0.4 0.4];
    case 'margin',      c = [0.7 0.7 0.7];
    case 'stub',        c = [0.7 0.7 0.7];
    otherwise
        if startsWith(cls, "tooth")
            k = sscanf(char(cls), 'tooth%d');
            c = tooth_colors(k, :);
        elseif startsWith(cls, "gap")
            c = [0.85 0.85 0.85];
        else
            c = [0 0 0];
        end
end
end

function p = ellipse_half_perim(a, b)
% Ramanujan's approximation for the half-perimeter of an ellipse.
p = pi * (3*(a+b) - sqrt((3*a+b)*(a+3*b))) / 2;
end

function b = solve_b_for_arc(L, a)
% Find b > 0 such that ellipse_half_perim(a, b) = L by bisection.
lo = 1e-6;
hi = max(1, 4*L/pi);
while ellipse_half_perim(a, hi) < L
    hi = 2*hi;
end
for it = 1:80
    mid = (lo + hi)/2;
    if ellipse_half_perim(a, mid) < L
        lo = mid;
    else
        hi = mid;
    end
end
b = (lo + hi)/2;
end
