% dolphin.m - Stylized dolphin silhouette through the Tutte + shape-opt
% pipeline. Streamlined body with a moderate dorsal fin (upward protrusion),
% a pectoral fin (downward protrusion), and a V-notched tail fluke at the
% rear. Easiest of the animal-geometry tests: moderate fin neck widths,
% no extreme aspect ratios in the protrusions.

rng(2);
porder = 4;
n_bnd  = 240;
n_int  = 400;

% --- Control points, CCW from snout ---------------------------------------
ctrl = [
    1.00,  0.00;    % snout tip
    0.92,  0.06;
    0.82,  0.12;
    0.70,  0.17;
    0.55,  0.19;
    0.42,  0.18;    % pre-dorsal
    0.38,  0.18;    % dorsal fin base front
    0.30,  0.32;    % dorsal fin peak
    0.20,  0.18;    % dorsal fin base back
    0.08,  0.15;
   -0.05,  0.11;
   -0.20,  0.08;
   -0.35,  0.04;    % peduncle upper
   -0.48,  0.10;    % upper fluke tip
   -0.48, -0.10;    % lower fluke tip (smooth trailing edge, no notch)
   -0.35, -0.04;    % peduncle lower
   -0.20, -0.08;
   -0.05, -0.11;
    0.08, -0.13;
    0.18, -0.15;    % pectoral fin base back
    0.26, -0.30;    % pectoral fin tip
    0.36, -0.16;    % pectoral fin base front
    0.50, -0.16;
    0.65, -0.15;
    0.80, -0.12;
    0.90, -0.07;
];
n_ctrl = size(ctrl, 1);

% Periodic spline to dense boundary
kw = 3;
ctrl_ext = [ctrl(end-kw+1:end, :); ctrl; ctrl(1:kw, :)];
t_ext    = (-kw : n_ctrl + kw - 1)' / n_ctrl;
t_dense  = (0:n_bnd-1)' / n_bnd;
xb = interp1(t_ext, ctrl_ext(:,1), t_dense, 'spline');
yb = interp1(t_ext, ctrl_ext(:,2), t_dense, 'spline');
bnd = [xb, yb];

% --- Interior rejection sampling ------------------------------------------
x_min = min(xb) - 0.02;  x_max = max(xb) + 0.02;
y_min = min(yb) - 0.02;  y_max = max(yb) + 0.02;
int_pts = zeros(0, 2);
while size(int_pts, 1) < n_int
    pt = [x_min + (x_max - x_min) * rand(), ...
          y_min + (y_max - y_min) * rand()];
    if inpolygon(pt(1), pt(2), xb, yb)
        int_pts(end+1, :) = pt;  %#ok<AGROW>
    end
end

run_animal_pipeline('dolphin', bnd, int_pts, porder);
