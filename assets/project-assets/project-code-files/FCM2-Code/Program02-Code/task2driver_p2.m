% task2driver_p2.m
%
% Task 2: Interpolation of f1(x) = (x - rho)^d
%
% Subroutines used:
%   mesh_uniform, mesh_cheby1, mesh_cheby2, leja_order
%   product_eval, bary1_coeffs, bary1_eval_dp
%   bary2_eval, newton_coeffs, newton_eval
%   cond_kappa, lebesgue_lambda, hilbert_H
%   f1
%
% Subtask 1: Interpolation setup on uniform, Cheby1, Cheby2 meshes (m=9)
% Subtask 2: Conditioning - kappa(x,n,y), kappa(x,n,1), Lambda_n, H_n
% Subtask 3: Accuracy and stability of single precision forms vs exact
%
% Output:
%   Table 1: Lambda_n and H_n for all mesh types
%   Table 2: Max error for all forms on all mesh types
%   Figure 1: kappa(x,n,1) Lebesgue function, uniform vs Cheby1 vs Cheby2
%   Figure 2: kappa(x,n,y) condition number, uniform vs Cheby1 vs Cheby2
%   Figure 3: Pointwise error all forms, uniform mesh (Higham-style)
%   Figure 4: Pointwise error all forms, Cheby1 mesh (Higham-style)
%   Figure 5: Pointwise error all forms, Cheby2 mesh
%   Figure 6: Stability bound vs observed error, uniform mesh

clear; clc; close all;  % Clear workspace, command window, and close all figures

% ----------------  settings  ----------------
rho   = 2.0;    % Root parameter for f1
d     = 9;      % Degree of f1 (m=9, so n+1=10 points)
n     = d;      % Number of interpolation points minus 1
a     = 0.5;    % Left endpoint (contains rho=2)
b     = 3.5;    % Right endpoint
n_pts = 1000;   % Number of query points for evaluation and plotting
u_sp  = single(eps('single'));  % Unit roundoff for single precision

x_eval = linspace(a, b, n_pts)';  % Dense query grid over [a, b]

fprintf('\n============================================================\n');  % Print separator
fprintf('TASK 2: f1(x) = (x - %.4g)^%d\n', rho, d);                          % Print task header
fprintf('Interval [%.4g, %.4g],  n=%d,  n_pts=%d\n', a, b, n, n_pts);        % Print settings
fprintf('============================================================\n\n');   % Print separator

% ============================================================
% SETUP: exact reference values via product form (double)
% ============================================================

f_exact = product_eval(double(x_eval), 1.0, rho*ones(d,1));  % Exact f1 at query points (double)

% ============================================================
% LOCAL HELPER: run all forms on a given mesh, return errors
% ============================================================

function [p_b2, p_inc, p_dec, p_lej] = eval_all_forms(x_nodes, rho, d, x_eval, mesh_type)
    % Evaluate all four single-precision interpolation forms on x_nodes.
    x_s   = single(x_nodes);                             % Cast nodes to single
    fv    = f1(x_s, rho, d);                             % Function values at nodes (single)
    beta  = bary2_coeffs(double(x_nodes), length(x_nodes)-1, mesh_type); % Bary2 weights O(n)

    p_b2  = bary2_eval(x_s, fv, single(beta), single(x_eval));  % Barycentric Form 2

    c_inc = newton_coeffs(x_s, fv);                              % Newton coeffs increasing
    p_inc = newton_eval(x_s, c_inc, single(x_eval));             % Newton increasing

    x_dec = flipud(x_s);  fv_dec = flipud(fv);                   % Reverse order
    c_dec = newton_coeffs(x_dec, fv_dec);                        % Newton coeffs decreasing
    p_dec = newton_eval(x_dec, c_dec, single(x_eval));           % Newton decreasing

    x_lej = single(leja_order(x_nodes));                         % Leja ordering
    fv_lej = f1(x_lej, rho, d);                                  % Function values at Leja nodes
    c_lej  = newton_coeffs(x_lej, fv_lej);                       % Newton coeffs Leja
    p_lej  = newton_eval(x_lej, c_lej, single(x_eval));          % Newton Leja
end

function [kappa_y, kappa_1, Lambda_n, H_n] = run_conditioning(x_nodes, rho, d, x_eval)
    % Compute condition numbers and summary statistics for a given mesh.
    fv_dp  = (double(x_nodes) - rho) .^ d;                     % Function values in double
    gamma  = bary1_coeffs(double(x_nodes), length(x_nodes)-1, ...
                          @(t) (t - rho).^d);                   % Bary1 weights (double)
    denom  = abs(product_eval(double(x_eval), 1.0, rho*ones(d,1)));  % Exact denominator

    [kappa_y, kappa_1] = cond_kappa(double(x_nodes), fv_dp, gamma, ...
                                     double(x_eval), denom);    % Condition numbers

    Lambda_n = lebesgue_lambda(kappa_1);  % Lebesgue constant
    H_n      = hilbert_H(kappa_y);        % H_n norm
end

% ============================================================
% SUBTASK 1 & 2: Meshes and conditioning
% ============================================================

fprintf('------------------------------------------------------------\n');  % Print separator
fprintf('SUBTASK 1 & 2: Mesh setup and conditioning\n');                    % Print section header
fprintf('------------------------------------------------------------\n\n'); % Print separator

x_uni  = mesh_uniform(a, b, n);   % Uniform mesh
x_ch1  = mesh_cheby1(a, b, n);    % Chebyshev first kind mesh
x_ch2  = mesh_cheby2(a, b, n);    % Chebyshev second kind mesh

[ky_uni, k1_uni, Lam_uni, H_uni] = run_conditioning(x_uni, rho, d, x_eval);   % Uniform conditioning
[ky_ch1, k1_ch1, Lam_ch1, H_ch1] = run_conditioning(x_ch1, rho, d, x_eval);   % Cheby1 conditioning
[ky_ch2, k1_ch2, Lam_ch2, H_ch2] = run_conditioning(x_ch2, rho, d, x_eval);   % Cheby2 conditioning

fprintf('Table 1: Conditioning summary for f1, n=%d\n', n);                                % Table header
fprintf('%-12s %-14s %-14s %-14s %-14s\n', 'Mesh', 'Lambda_n', 'H_n', ...                 % Column headers
    'max kappa_y', 'max kappa_1');
fprintf('%s\n', repmat('-',1,60));                                                           % Separator
fprintf('%-12s %-14.4e %-14.4e %-14.4e %-14.4e\n', 'Uniform', ...                         % Uniform row
    Lam_uni, H_uni, max(ky_uni(isfinite(ky_uni))), max(k1_uni));
fprintf('%-12s %-14.4e %-14.4e %-14.4e %-14.4e\n', 'Cheby1',  ...                         % Cheby1 row
    Lam_ch1, H_ch1, max(ky_ch1(isfinite(ky_ch1))), max(k1_ch1));
fprintf('%-12s %-14.4e %-14.4e %-14.4e %-14.4e\n', 'Cheby2',  ...                         % Cheby2 row
    Lam_ch2, H_ch2, max(ky_ch2(isfinite(ky_ch2))), max(k1_ch2));
fprintf('\n');  % Blank line

% Figure 1: Lebesgue function all meshes
figure(1); clf; hold on; grid on;                                                              % Create figure
plot(x_eval, k1_uni, 'LineWidth', 1.5, 'DisplayName', sprintf('Uniform (\\Lambda_n=%.2e)', Lam_uni));   % Uniform
plot(x_eval, k1_ch1, 'LineWidth', 1.5, 'DisplayName', sprintf('Cheby1  (\\Lambda_n=%.2e)', Lam_ch1));   % Cheby1
plot(x_eval, k1_ch2, 'LineWidth', 1.5, 'DisplayName', sprintf('Cheby2  (\\Lambda_n=%.2e)', Lam_ch2));   % Cheby2
xlabel('x'); ylabel('\kappa(x,n,1)');                                                          % Label axes
title(sprintf('Figure 1: Lebesgue function, f1, n=%d', n));                                   % Set title
legend('Location','best'); hold off;                                                           % Add legend

% Figure 2: kappa(x,n,y) all meshes
figure(2); clf; hold on; grid on;                                                              % Create figure
semilogy(x_eval, max(ky_uni, 1e-20), 'LineWidth', 1.5, 'DisplayName', sprintf('Uniform (H_n=%.2e)', H_uni));  % Uniform
semilogy(x_eval, max(ky_ch1, 1e-20), 'LineWidth', 1.5, 'DisplayName', sprintf('Cheby1  (H_n=%.2e)', H_ch1));  % Cheby1
semilogy(x_eval, max(ky_ch2, 1e-20), 'LineWidth', 1.5, 'DisplayName', sprintf('Cheby2  (H_n=%.2e)', H_ch2));  % Cheby2
xlabel('x'); ylabel('\kappa(x,n,y)');                                                          % Label axes
title(sprintf('Figure 2: Condition number kappa(x,n,y), f1, n=%d', n));                       % Set title
legend('Location','best'); hold off;                                                           % Add legend

% ============================================================
% SUBTASK 3: Accuracy and stability of single precision forms
% ============================================================

fprintf('------------------------------------------------------------\n');   % Print separator
fprintf('SUBTASK 3: Accuracy and stability of single precision forms\n');    % Print section header
fprintf('------------------------------------------------------------\n\n'); % Print separator

[p_b2_uni, p_inc_uni, p_dec_uni, p_lej_uni] = eval_all_forms(x_uni, rho, d, x_eval, 'uniform');   % Uniform forms
[p_b2_ch1, p_inc_ch1, p_dec_ch1, p_lej_ch1] = eval_all_forms(x_ch1, rho, d, x_eval, 'cheby1');   % Cheby1 forms
[p_b2_ch2, p_inc_ch2, p_dec_ch2, p_lej_ch2] = eval_all_forms(x_ch2, rho, d, x_eval, 'cheby2');   % Cheby2 forms

% Max errors for each form and mesh
err_b2_uni  = max(abs(double(p_b2_uni)  - f_exact));  % Bary2 uniform error
err_inc_uni = max(abs(double(p_inc_uni) - f_exact));  % Newton inc uniform error
err_dec_uni = max(abs(double(p_dec_uni) - f_exact));  % Newton dec uniform error
err_lej_uni = max(abs(double(p_lej_uni) - f_exact));  % Newton Leja uniform error

err_b2_ch1  = max(abs(double(p_b2_ch1)  - f_exact));  % Bary2 Cheby1 error
err_inc_ch1 = max(abs(double(p_inc_ch1) - f_exact));  % Newton inc Cheby1 error
err_dec_ch1 = max(abs(double(p_dec_ch1) - f_exact));  % Newton dec Cheby1 error
err_lej_ch1 = max(abs(double(p_lej_ch1) - f_exact));  % Newton Leja Cheby1 error

err_b2_ch2  = max(abs(double(p_b2_ch2)  - f_exact));  % Bary2 Cheby2 error
err_inc_ch2 = max(abs(double(p_inc_ch2) - f_exact));  % Newton inc Cheby2 error
err_dec_ch2 = max(abs(double(p_dec_ch2) - f_exact));  % Newton dec Cheby2 error
err_lej_ch2 = max(abs(double(p_lej_ch2) - f_exact));  % Newton Leja Cheby2 error

fprintf('Table 2: max |p(x) - f1(x)| for all forms and mesh types\n');                       % Table header
fprintf('%-24s %-14s %-14s %-14s\n', 'Form', 'Uniform', 'Cheby1', 'Cheby2');                 % Column headers
fprintf('%s\n', repmat('-',1,68));                                                             % Separator
fprintf('%-24s %-14.4e %-14.4e %-14.4e\n', 'Bary Form 2',       err_b2_uni,  err_b2_ch1,  err_b2_ch2);   % Bary2 row
fprintf('%-24s %-14.4e %-14.4e %-14.4e\n', 'Newton increasing', err_inc_uni, err_inc_ch1, err_inc_ch2);   % Newton inc row
fprintf('%-24s %-14.4e %-14.4e %-14.4e\n', 'Newton decreasing', err_dec_uni, err_dec_ch1, err_dec_ch2);   % Newton dec row
fprintf('%-24s %-14.4e %-14.4e %-14.4e\n', 'Newton Leja',       err_lej_uni, err_lej_ch1, err_lej_ch2);   % Newton Leja row
fprintf('\n');  % Blank line

% Stability bound for Bary Form 2 (Set 8 Section 5.4 / Higham 2004):
%   |p(x) - p_hat(x)| / |p(x)| <= (3n+4)*u*kappa_y + (3n+2)*u*kappa_1
% Multiply through by |f_exact| to get absolute bound on |error|.
u = double(u_sp);
stab_uni = ((3*n+4)*u*ky_uni + (3*n+2)*u*k1_uni) .* abs(f_exact);  % Bary2 stability bound uniform
stab_ch1 = ((3*n+4)*u*ky_ch1 + (3*n+2)*u*k1_ch1) .* abs(f_exact);  % Bary2 stability bound Cheby1
stab_ch2 = ((3*n+4)*u*ky_ch2 + (3*n+2)*u*k1_ch2) .* abs(f_exact);  % Bary2 stability bound Cheby2

% Figure 3: Higham-style error plot, uniform mesh
figure(3); clf; hold on; grid on;                                                              % Create figure
semilogy(x_eval, max(abs(double(p_b2_uni)  - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Bary2');        % Bary2
semilogy(x_eval, max(abs(double(p_inc_uni) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton inc');   % Newton inc
semilogy(x_eval, max(abs(double(p_dec_uni) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton dec');   % Newton dec
semilogy(x_eval, max(abs(double(p_lej_uni) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton Leja');  % Newton Leja
semilogy(x_eval, max(stab_uni, 1e-20), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Stability bound');               % Stability bound
xlabel('x'); ylabel('|p(x) - f1(x)|');                                                        % Label axes
title(sprintf('Figure 3: Pointwise error, f1, uniform mesh, n=%d', n));                       % Set title
legend('Location','best'); hold off;                                                           % Add legend

% Figure 4: Higham-style error plot, Cheby1 mesh
figure(4); clf; hold on; grid on;                                                              % Create figure
semilogy(x_eval, max(abs(double(p_b2_ch1)  - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Bary2');        % Bary2
semilogy(x_eval, max(abs(double(p_inc_ch1) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton inc');   % Newton inc
semilogy(x_eval, max(abs(double(p_dec_ch1) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton dec');   % Newton dec
semilogy(x_eval, max(abs(double(p_lej_ch1) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton Leja');  % Newton Leja
semilogy(x_eval, max(stab_ch1, 1e-20), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Stability bound');               % Stability bound
xlabel('x'); ylabel('|p(x) - f1(x)|');                                                        % Label axes
title(sprintf('Figure 4: Pointwise error, f1, Cheby1 mesh, n=%d', n));                        % Set title
legend('Location','best'); hold off;                                                           % Add legend

% Figure 5: error plot, Cheby2 mesh
figure(5); clf; hold on; grid on;                                                              % Create figure
semilogy(x_eval, max(abs(double(p_b2_ch2)  - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Bary2');        % Bary2
semilogy(x_eval, max(abs(double(p_inc_ch2) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton inc');   % Newton inc
semilogy(x_eval, max(abs(double(p_dec_ch2) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton dec');   % Newton dec
semilogy(x_eval, max(abs(double(p_lej_ch2) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton Leja');  % Newton Leja
semilogy(x_eval, max(stab_ch2, 1e-20), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Stability bound');               % Stability bound
xlabel('x'); ylabel('|p(x) - f1(x)|');                                                        % Label axes
title(sprintf('Figure 5: Pointwise error, f1, Cheby2 mesh, n=%d', n));                        % Set title
legend('Location','best'); hold off;                                                           % Add legend

% Figure 6: stability bound vs observed error on uniform mesh
figure(6); clf; hold on; grid on;                                                              % Create figure
semilogy(x_eval, max(abs(double(p_b2_uni) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Observed error (Bary2)');   % Observed error
semilogy(x_eval, max(stab_uni, 1e-20), 'k--', 'LineWidth', 1.5, 'DisplayName', '(3n+4)u\kappa_y + (3n+2)u\kappa_1 (Bary2 bound)');               % Stability bound
xlabel('x'); ylabel('error / bound');                                                          % Label axes
title(sprintf('Figure 6: Stability bound vs observed error, uniform, n=%d', n));              % Set title
legend('Location','best'); hold off;                                                           % Add legend

fprintf('DONE. Produced Tables 1-2 and Figures 1-6.\n');  % Print completion message
