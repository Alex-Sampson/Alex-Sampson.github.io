% task1driver_p3.m
%
% Task 1: Empirical validation of all P3 routines.
%
% Works through four validation experiments in order:
%   1. Polynomial reproduction -- splines and piecewise poly must reproduce
%      any cubic exactly (up to rounding) when given exact BCs.
%   2. Boundary condition checks -- both spline codes under both BC types
%      must give consistent results and match each other where expected.
%   3. Cross-code comparison -- spline1 vs spline2, and both vs piecewise
%      polynomial, on the same mesh and test functions.
%   4. Convergence study -- error vs number of subintervals M for all methods
%      on the Runge function.
%
% Subroutines exercised:
%   piecewise_interp, piecewise_eval
%   spline1_coeffs, spline1_eval
%   spline2_coeffs, spline2_eval
%   mesh_uniform
%
% Outputs:
%   Table 1:  Polynomial reproduction errors (clamped BCs, degrees 1-3)
%   Table 2:  BC consistency -- spline1 vs spline2 under both BC types
%   Table 3:  Cross-code max errors on Runge function for fixed M
%   Table 4:  Convergence table -- max error vs M for all methods
%   Figure 1: BC effect on g3 -- top: curves, bottom: errors
%   Figure 2: Pointwise error curves for all methods on Runge (fixed M)
%   Figure 3: Convergence plot -- all methods, log(error) vs log(h)
%   Figure 4: Convergence plot -- splines only, BC type comparison

clear; clc; close all;

% ---- color palette (clr_ prefix avoids clashes with output variable names) ----
clr_s1nat  = [0.122, 0.471, 0.706];   % steel blue   -- spline1 natural
clr_s1clmp = [0.651, 0.808, 0.890];   % light blue   -- spline1 clamped
clr_s2nat  = [0.890, 0.102, 0.110];   % crimson      -- spline2 natural
clr_s2clmp = [0.984, 0.604, 0.600];   % salmon       -- spline2 clamped
clr_pw1    = [0.200, 0.627, 0.173];   % forest green -- PW d=1
clr_pw2    = [0.596, 0.306, 0.639];   % purple       -- PW d=2
clr_pw3    = [1.000, 0.498, 0.000];   % orange       -- PW d=3
clr_herm   = [0.694, 0.349, 0.157];   % brown        -- PW Hermite
clr_exact  = [0.000, 0.000, 0.000];   % black        -- exact

% ---- global settings ----
a      = -1.0;
b      =  1.0;
n_pts  = 1000;
M_fix  = 8;

x_eval = linspace(a, b, n_pts)';

% function handles
f_runge  = @(x) 1 ./ (1 + 25*x.^2);
df_runge = @(x) -50*x ./ (1 + 25*x.^2).^2;

% test polynomials and their derivatives
g1  = @(x) 2*x + 1;
dg1 = @(x) 2 + 0*x;
g2  = @(x) x.^2 - x + 3;
dg2 = @(x) 2*x - 1;
g3  = @(x) x.^3 - 2*x.^2 + x - 1;
dg3 = @(x) 3*x.^2 - 4*x + 1;

% exact reference values on dense grid
f_exact  = f_runge(x_eval);
g1_exact = g1(x_eval);
g2_exact = g2(x_eval);
g3_exact = g3(x_eval);

fprintf('\n============================================================\n');
fprintf('TASK 1: Empirical Validation\n');
fprintf('Interval [%.4g, %.4g],  M_fix=%d subintervals,  %d query pts\n', a, b, M_fix, n_pts);
fprintf('============================================================\n\n');

% ============================================================
% VALIDATION 1: Polynomial reproduction
% ============================================================

fprintf('------------------------------------------------------------\n');
fprintf('VALIDATION 1: Polynomial reproduction (clamped BCs)\n');
fprintf('------------------------------------------------------------\n\n');

x_knots = mesh_uniform(a, b, M_fix);
f_g1 = g1(x_knots);
f_g2 = g2(x_knots);
f_g3 = g3(x_knots);

da1 = dg1(a);  db1 = dg1(b);
da2 = dg2(a);  db2 = dg2(b);
da3 = dg3(a);  db3 = dg3(b);

% spline1 clamped
[s2_g1, br1, fv1] = spline1_coeffs(x_knots, f_g1, 'clamped', [da1, db1]);
[s2_g2, br2, fv2] = spline1_coeffs(x_knots, f_g2, 'clamped', [da2, db2]);
[s2_g3, br3, fv3] = spline1_coeffs(x_knots, f_g3, 'clamped', [da3, db3]);

p_s1_g1 = zeros(n_pts,1); p_s1_g2 = zeros(n_pts,1); p_s1_g3 = zeros(n_pts,1);
for k = 1 : n_pts
    p_s1_g1(k) = spline1_eval(s2_g1, br1, fv1, x_eval(k));
    p_s1_g2(k) = spline1_eval(s2_g2, br2, fv2, x_eval(k));
    p_s1_g3(k) = spline1_eval(s2_g3, br3, fv3, x_eval(k));
end

% spline2 clamped
[al_g1, br_b1] = spline2_coeffs(x_knots, f_g1, 'clamped', [da1, db1]);
[al_g2, br_b2] = spline2_coeffs(x_knots, f_g2, 'clamped', [da2, db2]);
[al_g3, br_b3] = spline2_coeffs(x_knots, f_g3, 'clamped', [da3, db3]);

p_s2_g1 = zeros(n_pts,1); p_s2_g2 = zeros(n_pts,1); p_s2_g3 = zeros(n_pts,1);
for k = 1 : n_pts
    p_s2_g1(k) = spline2_eval(al_g1, br_b1, x_eval(k));
    p_s2_g2(k) = spline2_eval(al_g2, br_b2, x_eval(k));
    p_s2_g3(k) = spline2_eval(al_g3, br_b3, x_eval(k));
end

% piecewise poly (coef_ prefix avoids color variable name collision)
[coef_pw1, bk_pw1, nd_pw1] = piecewise_interp(g1, [], a, b, M_fix, 1, 'uniform', false);
[coef_pw2, bk_pw2, nd_pw2] = piecewise_interp(g2, [], a, b, M_fix, 2, 'uniform', false);
[coef_pw3, bk_pw3, nd_pw3] = piecewise_interp(g3, [], a, b, M_fix, 3, 'uniform', false);

p_pw1 = zeros(n_pts,1); p_pw2 = zeros(n_pts,1); p_pw3 = zeros(n_pts,1);
for k = 1 : n_pts
    p_pw1(k) = piecewise_eval(coef_pw1, bk_pw1, nd_pw1, x_eval(k), false);
    p_pw2(k) = piecewise_eval(coef_pw2, bk_pw2, nd_pw2, x_eval(k), false);
    p_pw3(k) = piecewise_eval(coef_pw3, bk_pw3, nd_pw3, x_eval(k), false);
end

err_s1 = [max(abs(p_s1_g1-g1_exact)), max(abs(p_s1_g2-g2_exact)), max(abs(p_s1_g3-g3_exact))];
err_s2 = [max(abs(p_s2_g1-g1_exact)), max(abs(p_s2_g2-g2_exact)), max(abs(p_s2_g3-g3_exact))];
err_pw = [max(abs(p_pw1-g1_exact)),   max(abs(p_pw2-g2_exact)),   max(abs(p_pw3-g3_exact))];

fprintf('Table 1: Polynomial reproduction -- max |approx - exact| (clamped BC, M=%d)\n', M_fix);
fprintf('%-20s %-18s %-18s %-18s\n', 'Method', 'g1 (deg 1)', 'g2 (deg 2)', 'g3 (deg 3)');
fprintf('%s\n', repmat('-',1,76));
fprintf('%-20s %-18.4e %-18.4e %-18.4e\n', 'Spline 1',        err_s1(1), err_s1(2), err_s1(3));
fprintf('%-20s %-18.4e %-18.4e %-18.4e\n', 'Spline 2',        err_s2(1), err_s2(2), err_s2(3));
fprintf('%-20s %-18.4e %-18.4e %-18.4e\n', 'Piecewise d=deg', err_pw(1), err_pw(2), err_pw(3));
fprintf('\n');

assert(err_s1(1) < 1e-10, 'Spline1 failed to reproduce linear g1');
assert(err_s1(2) < 1e-10, 'Spline1 failed to reproduce quadratic g2');
assert(err_s1(3) < 1e-10, 'Spline1 failed to reproduce cubic g3');
assert(err_s2(1) < 1e-10, 'Spline2 failed to reproduce linear g1');
assert(err_s2(2) < 1e-10, 'Spline2 failed to reproduce quadratic g2');
assert(err_s2(3) < 1e-10, 'Spline2 failed to reproduce cubic g3');
assert(err_pw(1) < 1e-10, 'Piecewise d=1 failed to reproduce linear g1');
assert(err_pw(2) < 1e-10, 'Piecewise d=2 failed to reproduce quadratic g2');
assert(err_pw(3) < 1e-10, 'Piecewise d=3 failed to reproduce cubic g3');
fprintf('All polynomial reproduction assertions passed.\n\n');

% ============================================================
% VALIDATION 2: Boundary condition consistency
% ============================================================

fprintf('------------------------------------------------------------\n');
fprintf('VALIDATION 2: Boundary condition consistency\n');
fprintf('------------------------------------------------------------\n\n');

[s2_nat,  br_n, fv_n] = spline1_coeffs(x_knots, f_g3, 'natural', []);
[s2_clmp, br_c, fv_c] = spline1_coeffs(x_knots, f_g3, 'clamped', [da3, db3]);
[al_nat,  br_bn]      = spline2_coeffs(x_knots, f_g3, 'natural', []);
[al_clmp, br_bc]      = spline2_coeffs(x_knots, f_g3, 'clamped', [da3, db3]);

p_s1_nat  = zeros(n_pts,1); p_s1_clmp = zeros(n_pts,1);
p_s2_nat  = zeros(n_pts,1); p_s2_clmp = zeros(n_pts,1);
for k = 1 : n_pts
    p_s1_nat(k)  = spline1_eval(s2_nat,  br_n,  fv_n,  x_eval(k));
    p_s1_clmp(k) = spline1_eval(s2_clmp, br_c,  fv_c,  x_eval(k));
    p_s2_nat(k)  = spline2_eval(al_nat,  br_bn, x_eval(k));
    p_s2_clmp(k) = spline2_eval(al_clmp, br_bc, x_eval(k));
end

e_s1_nat        = max(abs(p_s1_nat  - g3_exact));
e_s1_clmp       = max(abs(p_s1_clmp - g3_exact));
e_s2_nat        = max(abs(p_s2_nat  - g3_exact));
e_s2_clmp       = max(abs(p_s2_clmp - g3_exact));
e_s1_vs_s2_nat  = max(abs(p_s1_nat  - p_s2_nat));
e_s1_vs_s2_clmp = max(abs(p_s1_clmp - p_s2_clmp));

fprintf('Table 2: BC consistency on cubic g3\n');
fprintf('%-28s %-18s %-18s\n', 'Comparison', 'Natural BC', 'Clamped BC');
fprintf('%s\n', repmat('-',1,66));
fprintf('%-28s %-18.4e %-18.4e\n', 'Spline1 vs exact',   e_s1_nat,        e_s1_clmp);
fprintf('%-28s %-18.4e %-18.4e\n', 'Spline2 vs exact',   e_s2_nat,        e_s2_clmp);
fprintf('%-28s %-18.4e %-18.4e\n', 'Spline1 vs Spline2', e_s1_vs_s2_nat,  e_s1_vs_s2_clmp);
fprintf('\n');

assert(e_s1_clmp       < 1e-10, 'Spline1 clamped failed on cubic g3');
assert(e_s2_clmp       < 1e-10, 'Spline2 clamped failed on cubic g3');
assert(e_s1_vs_s2_nat  < 1e-10, 'Spline1 and Spline2 disagree under natural BC');
assert(e_s1_vs_s2_clmp < 1e-10, 'Spline1 and Spline2 disagree under clamped BC');
fprintf('All BC consistency assertions passed.\n\n');

% ---- Figure 1: BC effect on g3 ----
% top panel: overlay all four splines on exact g3 so endpoint deviation
% from natural BCs is visible in the function itself.
% bottom panel: log-scale errors make it unambiguous that clamped sits
% at machine precision while natural has a real O(1) error near endpoints.
figure(1); clf;

subplot(2,1,1);
hold on; grid on; box on;
plot(x_eval,  g3_exact,   '-',  'Color', clr_exact,  'LineWidth', 2.0, 'DisplayName', 'Exact g_3');
plot(x_eval,  p_s1_clmp,  '--', 'Color', clr_s1nat,  'LineWidth', 1.8, 'DisplayName', 'S1 clamped');
plot(x_eval,  p_s1_nat,   ':',  'Color', clr_s1nat,  'LineWidth', 1.8, 'DisplayName', 'S1 natural');
plot(x_eval,  p_s2_clmp,  '--', 'Color', clr_s2nat,  'LineWidth', 1.8, 'DisplayName', 'S2 clamped');
plot(x_eval,  p_s2_nat,   ':',  'Color', clr_s2nat,  'LineWidth', 1.8, 'DisplayName', 'S2 natural');
plot(x_knots, f_g3, 'ko', 'MarkerSize', 5, 'HandleVisibility', 'off');
xlabel('x'); ylabel('s(x)');
title('Figure 1 (top): Spline fits on g_3(x) -- clamped (dashed) vs natural (dotted)');
legend('Location', 'best'); hold off;

subplot(2,1,2);
hold on; grid on; box on;
semilogy(x_eval, max(abs(p_s1_clmp - g3_exact), 1e-18), '--', 'Color', clr_s1nat,  'LineWidth', 1.8, 'DisplayName', 'S1 clamped');
semilogy(x_eval, max(abs(p_s1_nat  - g3_exact), 1e-18), ':',  'Color', clr_s1nat,  'LineWidth', 1.8, 'DisplayName', 'S1 natural');
semilogy(x_eval, max(abs(p_s2_clmp - g3_exact), 1e-18), '--', 'Color', clr_s2nat,  'LineWidth', 1.8, 'DisplayName', 'S2 clamped');
semilogy(x_eval, max(abs(p_s2_nat  - g3_exact), 1e-18), ':',  'Color', clr_s2nat,  'LineWidth', 1.8, 'DisplayName', 'S2 natural');
xlabel('x'); ylabel('|s(x) - g_3(x)|');
title('Figure 1 (bottom): Clamped BCs reproduce g_3 to machine precision; natural BCs do not');
legend('Location', 'best'); hold off;

% ============================================================
% VALIDATION 3: Cross-code comparison on Runge function
% ============================================================

fprintf('------------------------------------------------------------\n');
fprintf('VALIDATION 3: Cross-code comparison on Runge function (M=%d)\n', M_fix);
fprintf('------------------------------------------------------------\n\n');

x_knots_r = mesh_uniform(a, b, M_fix);
f_knots_r  = f_runge(x_knots_r);
da_r = df_runge(a);
db_r = df_runge(b);

[s2_rn, br_rn, fv_rn] = spline1_coeffs(x_knots_r, f_knots_r, 'natural', []);
[s2_rc, br_rc, fv_rc] = spline1_coeffs(x_knots_r, f_knots_r, 'clamped', [da_r, db_r]);
[al_rn, br_rbn]       = spline2_coeffs(x_knots_r, f_knots_r, 'natural', []);
[al_rc, br_rbc]       = spline2_coeffs(x_knots_r, f_knots_r, 'clamped', [da_r, db_r]);

[coef_p1, bk_p1, nd_p1] = piecewise_interp(f_runge, [],       a, b, M_fix, 1, 'uniform', false);
[coef_p2, bk_p2, nd_p2] = piecewise_interp(f_runge, [],       a, b, M_fix, 2, 'uniform', false);
[coef_p3, bk_p3, nd_p3] = piecewise_interp(f_runge, [],       a, b, M_fix, 3, 'uniform', false);
[coef_ph, bk_ph, nd_ph] = piecewise_interp(f_runge, df_runge, a, b, M_fix, 3, 'uniform', true);

p_s1n  = zeros(n_pts,1); p_s1c  = zeros(n_pts,1);
p_s2n  = zeros(n_pts,1); p_s2c  = zeros(n_pts,1);
p_d1   = zeros(n_pts,1); p_d2   = zeros(n_pts,1);
p_d3   = zeros(n_pts,1); p_herm = zeros(n_pts,1);

for k = 1 : n_pts
    p_s1n(k)  = spline1_eval(s2_rn,  br_rn,  fv_rn,  x_eval(k));
    p_s1c(k)  = spline1_eval(s2_rc,  br_rc,  fv_rc,  x_eval(k));
    p_s2n(k)  = spline2_eval(al_rn,  br_rbn, x_eval(k));
    p_s2c(k)  = spline2_eval(al_rc,  br_rbc, x_eval(k));
    p_d1(k)   = piecewise_eval(coef_p1, bk_p1, nd_p1, x_eval(k), false);
    p_d2(k)   = piecewise_eval(coef_p2, bk_p2, nd_p2, x_eval(k), false);
    p_d3(k)   = piecewise_eval(coef_p3, bk_p3, nd_p3, x_eval(k), false);
    p_herm(k) = piecewise_eval(coef_ph, bk_ph, nd_ph, x_eval(k), true);
end

methods  = {'S1 natural','S1 clamped','S2 natural','S2 clamped', ...
            'PW d=1','PW d=2','PW d=3','PW Hermite'};
all_p    = {p_s1n, p_s1c, p_s2n, p_s2c, p_d1, p_d2, p_d3, p_herm};
all_errs = cellfun(@(p) max(abs(p - f_exact)), all_p);

fprintf('Table 3: Max |approx - f_runge| on [%.4g,%.4g], M=%d, uniform knots\n', a, b, M_fix);
fprintf('%-18s %-18s\n', 'Method', 'max error');
fprintf('%s\n', repmat('-',1,38));
for m = 1 : length(methods)
    fprintf('%-18s %-18.4e\n', methods{m}, all_errs(m));
end
fprintf('\n');

% ---- Figure 2: pointwise errors on Runge ----
% semilogy reveals order-of-magnitude differences between methods clearly.
% blue family = S1, red family = S2, remaining colors = PW variants.
% solid = natural BC, dashed = clamped BC within each spline family.
figure(2); clf; hold on; grid on; box on;
semilogy(x_eval, max(abs(p_s1n  - f_exact), 1e-20), '-',  'Color', clr_s1nat,  'LineWidth', 1.8, 'DisplayName', 'S1 natural');
semilogy(x_eval, max(abs(p_s1c  - f_exact), 1e-20), '--', 'Color', clr_s1nat,  'LineWidth', 1.8, 'DisplayName', 'S1 clamped');
semilogy(x_eval, max(abs(p_s2n  - f_exact), 1e-20), '-',  'Color', clr_s2nat,  'LineWidth', 1.8, 'DisplayName', 'S2 natural');
semilogy(x_eval, max(abs(p_s2c  - f_exact), 1e-20), '--', 'Color', clr_s2nat,  'LineWidth', 1.8, 'DisplayName', 'S2 clamped');
semilogy(x_eval, max(abs(p_d1   - f_exact), 1e-20), '-',  'Color', clr_pw1,    'LineWidth', 1.8, 'DisplayName', 'PW d=1');
semilogy(x_eval, max(abs(p_d2   - f_exact), 1e-20), '-',  'Color', clr_pw2,    'LineWidth', 1.8, 'DisplayName', 'PW d=2');
semilogy(x_eval, max(abs(p_d3   - f_exact), 1e-20), '-',  'Color', clr_pw3,    'LineWidth', 1.8, 'DisplayName', 'PW d=3');
semilogy(x_eval, max(abs(p_herm - f_exact), 1e-20), '-',  'Color', clr_herm,   'LineWidth', 1.8, 'DisplayName', 'PW Hermite');
xlabel('x'); ylabel('|approx(x) - f(x)|');
title(sprintf('Figure 2: Pointwise errors on Runge function, M=%d uniform subintervals', M_fix));
legend('Location', 'best'); hold off;

% ============================================================
% VALIDATION 4: Convergence study
% ============================================================

fprintf('------------------------------------------------------------\n');
fprintf('VALIDATION 4: Convergence study vs number of subintervals M\n');
fprintf('------------------------------------------------------------\n\n');

M_vals   = [2, 4, 8, 16, 32, 64];
nM       = length(M_vals);
err_conv = zeros(8, nM);

for im = 1 : nM
    M    = M_vals(im);
    xk   = mesh_uniform(a, b, M);
    fk   = f_runge(xk);
    da_c = df_runge(a);
    db_c = df_runge(b);

    [s2n, brn, fvn]   = spline1_coeffs(xk, fk, 'natural', []);
    [s2c, brc, fvc]   = spline1_coeffs(xk, fk, 'clamped', [da_c, db_c]);
    [aln, bran]       = spline2_coeffs(xk, fk, 'natural', []);
    [alc, brac]       = spline2_coeffs(xk, fk, 'clamped', [da_c, db_c]);
    [cp1, bkp1, ndp1] = piecewise_interp(f_runge, [],       a, b, M, 1, 'uniform', false);
    [cp2, bkp2, ndp2] = piecewise_interp(f_runge, [],       a, b, M, 2, 'uniform', false);
    [cp3, bkp3, ndp3] = piecewise_interp(f_runge, [],       a, b, M, 3, 'uniform', false);
    [cph, bkph, ndph] = piecewise_interp(f_runge, df_runge, a, b, M, 3, 'uniform', true);

    ps1n = zeros(n_pts,1); ps1c = zeros(n_pts,1);
    ps2n = zeros(n_pts,1); ps2c = zeros(n_pts,1);
    pd1  = zeros(n_pts,1); pd2  = zeros(n_pts,1);
    pd3  = zeros(n_pts,1); pdh  = zeros(n_pts,1);

    for k = 1 : n_pts
        ps1n(k) = spline1_eval(s2n,  brn,  fvn,  x_eval(k));
        ps1c(k) = spline1_eval(s2c,  brc,  fvc,  x_eval(k));
        ps2n(k) = spline2_eval(aln,  bran, x_eval(k));
        ps2c(k) = spline2_eval(alc,  brac, x_eval(k));
        pd1(k)  = piecewise_eval(cp1, bkp1, ndp1, x_eval(k), false);
        pd2(k)  = piecewise_eval(cp2, bkp2, ndp2, x_eval(k), false);
        pd3(k)  = piecewise_eval(cp3, bkp3, ndp3, x_eval(k), false);
        pdh(k)  = piecewise_eval(cph, bkph, ndph, x_eval(k), true);
    end

    err_conv(1,im) = max(abs(ps1n - f_exact));
    err_conv(2,im) = max(abs(ps1c - f_exact));
    err_conv(3,im) = max(abs(ps2n - f_exact));
    err_conv(4,im) = max(abs(ps2c - f_exact));
    err_conv(5,im) = max(abs(pd1  - f_exact));
    err_conv(6,im) = max(abs(pd2  - f_exact));
    err_conv(7,im) = max(abs(pd3  - f_exact));
    err_conv(8,im) = max(abs(pdh  - f_exact));
end

conv_names = {'S1 natural','S1 clamped','S2 natural','S2 clamped', ...
              'PW d=1','PW d=2','PW d=3','PW Hermite'};

fprintf('Table 4: Max error vs M (Runge function, uniform mesh)\n');
header = sprintf('%-14s', 'Method');
for im = 1 : nM
    header = [header, sprintf('%-12s', ['M=', num2str(M_vals(im))])]; %#ok<AGROW>
end
fprintf('%s\n', header);
fprintf('%s\n', repmat('-', 1, 14 + 12*nM));
for m = 1 : 8
    row = sprintf('%-14s', conv_names{m});
    for im = 1 : nM
        row = [row, sprintf('%-12.3e', err_conv(m,im))]; %#ok<AGROW>
    end
    fprintf('%s\n', row);
end
fprintf('\n');

h_vals = (b - a) ./ M_vals;
h_ref  = linspace(min(h_vals), max(h_vals), 50);

% ---- Figure 3: all methods convergence log-log ----
% every method on one plot with reference slope lines.
% log-log is the standard way to read off O(h^p) rates visually --
% each method appears as a straight line whose slope is p.
figure(3); clf; hold on; grid on; box on;
loglog(h_vals, err_conv(1,:), '-o',  'Color', clr_s1nat,  'LineWidth', 1.8, 'DisplayName', 'S1 natural');
loglog(h_vals, err_conv(2,:), '--o', 'Color', clr_s1nat,  'LineWidth', 1.8, 'DisplayName', 'S1 clamped');
loglog(h_vals, err_conv(3,:), '-s',  'Color', clr_s2nat,  'LineWidth', 1.8, 'DisplayName', 'S2 natural');
loglog(h_vals, err_conv(4,:), '--s', 'Color', clr_s2nat,  'LineWidth', 1.8, 'DisplayName', 'S2 clamped');
loglog(h_vals, err_conv(5,:), '-^',  'Color', clr_pw1,    'LineWidth', 1.8, 'DisplayName', 'PW d=1');
loglog(h_vals, err_conv(6,:), '-^',  'Color', clr_pw2,    'LineWidth', 1.8, 'DisplayName', 'PW d=2');
loglog(h_vals, err_conv(7,:), '-^',  'Color', clr_pw3,    'LineWidth', 1.8, 'DisplayName', 'PW d=3');
loglog(h_vals, err_conv(8,:), '-d',  'Color', clr_herm,   'LineWidth', 1.8, 'DisplayName', 'PW Hermite');
loglog(h_ref,  0.3*h_ref.^2, 'k:',  'LineWidth', 1.2, 'DisplayName', 'O(h^2) ref');
loglog(h_ref,  0.3*h_ref.^3, 'k--', 'LineWidth', 1.2, 'DisplayName', 'O(h^3) ref');
loglog(h_ref,  0.3*h_ref.^4, 'k-',  'LineWidth', 1.2, 'DisplayName', 'O(h^4) ref');
xlabel('h = (b-a)/M  (subinterval width)'); ylabel('max |approx - f|');
title('Figure 3: Convergence rates -- all methods on Runge function');
legend('Location', 'southeast'); hold off;

% ---- Figure 4: splines only -- BC type effect on convergence ----
% strips away the piecewise poly lines so the BC question is the sole focus.
% the key observation: clamped BCs achieve O(h^4) while natural BCs
% typically achieve only O(h^3) on smooth functions where s''(endpoints) != 0.
% S1 and S2 should track each other exactly for each BC type.
figure(4); clf; hold on; grid on; box on;
loglog(h_vals, err_conv(1,:), '-o',  'Color', clr_s1nat,  'LineWidth', 1.8, 'DisplayName', 'S1 natural');
loglog(h_vals, err_conv(2,:), '--o', 'Color', clr_s1nat,  'LineWidth', 1.8, 'DisplayName', 'S1 clamped');
loglog(h_vals, err_conv(3,:), '-s',  'Color', clr_s2nat,  'LineWidth', 1.8, 'DisplayName', 'S2 natural');
loglog(h_vals, err_conv(4,:), '--s', 'Color', clr_s2nat,  'LineWidth', 1.8, 'DisplayName', 'S2 clamped');
loglog(h_ref,  0.3*h_ref.^3, 'k--', 'LineWidth', 1.2, 'DisplayName', 'O(h^3) ref');
loglog(h_ref,  0.3*h_ref.^4, 'k-',  'LineWidth', 1.2, 'DisplayName', 'O(h^4) ref');
xlabel('h = (b-a)/M  (subinterval width)'); ylabel('max |s(x) - f(x)|');
title('Figure 4: Spline convergence -- clamped BCs give O(h^4), natural give O(h^3)');
legend('Location', 'southeast'); hold off;

fprintf('DONE. Produced Tables 1-4 and Figures 1-4.\n');