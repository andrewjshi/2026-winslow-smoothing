% w_regression.m - Regression test for the W-modified smoother.
%
% Verifies that elliptic_smoothing_W reproduces elliptic_smoothing
% exactly in two scenarios:
%   (a) msh.W_K absent          -> old code path, bit-for-bit match.
%   (b) msh.W_K = identity (per element) -> new code path with W = I,
%                                            should match up to roundoff.
% Uses the standard w2.m setup (unskewed 11x11 mesh, side-4 perturbation).

porder = 4;
n = 10;
msh = mshsquare(n+1, n+1);
msh = nodealloc(msh, porder);

[~, edg] = getbndnodes(msh, dginit(msh), 4);
x = msh.p1(:,1,:);
y = msh.p1(:,2,:);
newx = x;
newy = y;
newx(edg) = x(edg) + 0.5*sin(pi*y(edg));
curvep1 = cat(2, newx, newy);

doplot = false;

% --- Baseline: original smoother ---------------------------------------
msh_baseline = elliptic_smoothing(msh, curvep1, doplot);

% --- Scenario (a): modified smoother without W_h ------------------------
msh_mod_nowh = elliptic_smoothing_W(msh, curvep1, doplot);
err_a = norm(msh_baseline.p1(:) - msh_mod_nowh.p1(:), inf);
fprintf('\n[Scenario a] ||p1_base - p1_W(noWh)||_inf = %.3e\n', err_a);

% --- Scenario (b): modified smoother with W_h = I per node -------------
nt = size(msh.t, 2);
ns = size(msh.s, 1);
msh_W = msh;
msh_W.W_h = repmat(reshape(eye(2), 1, 2, 2), [ns, 1, 1, nt]);
msh_mod_I = elliptic_smoothing_W(msh_W, curvep1, doplot);
err_b = norm(msh_baseline.p1(:) - msh_mod_I.p1(:), inf);
fprintf('[Scenario b] ||p1_base - p1_W(Wh=I)||_inf  = %.3e\n', err_b);

% --- Verdict ------------------------------------------------------------
tol = 1e-12;
pass_a = err_a < tol;
pass_b = err_b < tol;
fprintf('\nScenario (a) %s  (tol %.0e)\n', pass_str(pass_a), tol);
fprintf('Scenario (b) %s  (tol %.0e)\n', pass_str(pass_b), tol);
if pass_a && pass_b
    fprintf('\nALL REGRESSION TESTS PASSED.\n');
else
    error('Regression tests failed.');
end


function s = pass_str(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
