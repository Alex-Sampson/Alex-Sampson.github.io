% task5driver_p2.m
%
% Task 5: Interpolation of f4(x) = 1 / (1 + 25*x^2)  (Runge's function)
%
% Key distinction from Tasks 2-4: f4 is not a polynomial, so no product
% form is available.  The "exact" reference and conditioning denominator
% are both computed via bary1_eval_dp in double precision throughout.
%
% Subroutines used:
%   mesh_uniform, mesh_cheby1, mesh_cheby2, leja_order
%   bary1_coeffs, bary1_eval_dp
%   bary2_eval, newton_coeffs, newton_eval
%   cond_kappa, lebesgue_lambda, hilbert_H
%   f4
%
% Subtask 1: Conditioning sweep over n for uniform, Cheby1, Cheby2 meshes
% Subtask 2: Accuracy and stability of single precision forms at n=29
%            Higham-style plots for uniform and Cheby1 meshes
% Subtask 3: Convergence investigation - ||f4 - p_n||_inf vs n
%            Uniform (Runge divergence), Cheby1 and Cheby2 (convergence)
%            Empirically determine n needed for target accuracy levels
%            Identify single precision floor
%
% Output:
%   Table 1: Lambda_n and H_n for all n and mesh types
%   Table 2: Max error for all forms, n=29, uniform and Cheby1
%   Table 3: Convergence table - inf-norm error vs n for Cheby1 and Cheby2
%   Figure 1: Lambda_n vs n for uniform, Cheby1, Cheby2
%   Figure 2: H_n vs n for uniform, Cheby1, Cheby2
%   Figure 3: Pointwise error all forms, uniform mesh, n=29 (Higham-style)
%   Figure 4: Pointwise error all forms, Cheby1 mesh, n=29 (Higham-style)
%   Figure 5: ||f4 - p_n||_inf vs n (convergence investigation)
%   Figure 6: Convergence zoom - Cheby1 and Cheby2 only, with precision floor

clear; clc; close all;  % Clear workspace, command window, and close all figures

% ----------------  settings  ----------------
a      = -1.0;   % Left endpoint (standard Runge interval)
b      =  1.0;   % Right endpoint
n_pts  = 1000;   % Number of query points for evaluation and plotting
u_sp   = single(eps('single'));  % Unit roundoff for single precision

n_sweep  = [5, 10, 15, 20, 25, 29, 35, 40, 50];  % Degrees for conditioning and convergence sweep
n_higham = 29;                                     % Degree for Higham-style plots (30 points)

% Target accuracy levels for convergence investigation (Subtask 3)
target_levels = [1e-2, 1e-4, 1e-6, 1e-8];  % Target inf-norm errors

x_eval = linspace(a, b, n_pts)';  % Dense query grid over [a, b]

% Exact f4 values at query points (double precision reference)
f_exact = double(f4(double(x_eval)));  % f4 in double precision serves as exact reference

fprintf('\n============================================================\n');  % Print separator
fprintf('TASK 5: f4(x) = 1 / (1 + 25*x^2)  (Runge function)\n');             % Print task header
fprintf('Interval [%.4g, %.4g]\n', a, b);                                     % Print interval
fprintf('n sweep: ');  fprintf('%d ', n_sweep);  fprintf('\n');                % Print degree list
fprintf('Note: no product form; bary1_eval_dp used as exact reference\n');    % Print key note
fprintf('============================================================\n\n');  % Print separator

% ============================================================
% LOCAL HELPERS
% ============================================================

function [p_b2, p_inc, p_dec, p_lej] = eval_all_forms(x_nodes, x_eval, mesh_type)
    % Evaluate all four single-precision forms for f4 on a given mesh.
    x_s   = single(x_nodes);                              % Cast nodes to single
    fv    = f4(x_s);                                      % f4 values at nodes (single)
    beta  = bary2_coeffs(double(x_nodes), length(x_nodes)-1, mesh_type); % Bary2 weights O(n)

    p_b2  = bary2_eval(x_s, fv, single(beta), single(x_eval));  % Barycentric Form 2

    c_inc = newton_coeffs(x_s, fv);                               % Newton coeffs increasing
    p_inc = newton_eval(x_s, c_inc, single(x_eval));              % Newton increasing

    x_dec = flipud(x_s);  fv_dec = flipud(fv);                    % Reverse order
    c_dec = newton_coeffs(x_dec, fv_dec);                         % Newton coeffs decreasing
    p_dec = newton_eval(x_dec, c_dec, single(x_eval));            % Newton decreasing

    x_lej = single(leja_order(x_nodes));                          % Leja ordering
    fv_lej = f4(x_lej);                                           % f4 values at Leja nodes
    c_lej  = newton_coeffs(x_lej, fv_lej);                        % Newton coeffs Leja
    p_lej  = newton_eval(x_lej, c_lej, single(x_eval));           % Newton Leja
end

function [kappa_y, kappa_1, Lambda_n, H_n, f_ref] = run_conditioning(x_nodes, x_eval)
    % Compute condition numbers and exact reference for f4 on a given mesh.
    % bary1_eval_dp provides both p_exact (f_ref) and the denominator.
    n      = length(x_nodes) - 1;                                  % Degree
    fv_dp  = double(f4(double(x_nodes)));                          % f4 values in double
    gamma  = bary1_coeffs(double(x_nodes), n, ...
                          @(t) double(f4(double(t))));              % Bary1 weights (double)

    [f_ref, ~, denom] = bary1_eval_dp(double(x_nodes), fv_dp, ...
                                       gamma, double(x_eval));      % Exact ref and denominator

    [kappa_y, kappa_1] = cond_kappa(double(x_nodes), fv_dp, gamma, ...
                                     double(x_eval), abs(denom));   % Condition numbers

    Lambda_n = lebesgue_lambda(kappa_1);  % Lebesgue constant
    H_n      = hilbert_H(kappa_y);        % H_n norm
end

% ============================================================
% SUBTASK 1: Conditioning sweep over all n
% ============================================================

fprintf('------------------------------------------------------------\n');   % Print separator
fprintf('SUBTASK 1: Conditioning sweep\n');                                   % Section header
fprintf('------------------------------------------------------------\n\n'); % Print separator

fprintf('Table 1: Conditioning summary for f4 (all n and meshes)\n');                           % Table header
fprintf('%-6s %-10s %-14s %-14s %-14s %-14s\n', ...                                             % Column headers
    'n', 'Mesh', 'Lambda_n', 'H_n', 'max kappa_y', 'max kappa_1');
fprintf('%s\n', repmat('-', 1, 76));  % Separator

Lam_uni_all = zeros(length(n_sweep), 1);  % Lambda_n storage uniform
Lam_ch1_all = zeros(length(n_sweep), 1);  % Lambda_n storage Cheby1
Lam_ch2_all = zeros(length(n_sweep), 1);  % Lambda_n storage Cheby2
H_uni_all   = zeros(length(n_sweep), 1);  % H_n storage uniform
H_ch1_all   = zeros(length(n_sweep), 1);  % H_n storage Cheby1
H_ch2_all   = zeros(length(n_sweep), 1);  % H_n storage Cheby2

for ni = 1 : length(n_sweep)  % Loop over each degree in sweep
    n = n_sweep(ni);           % Get current degree

    x_uni = mesh_uniform(a, b, n);  % Uniform mesh
    x_ch1 = mesh_cheby1(a, b, n);  % Cheby1 mesh
    x_ch2 = mesh_cheby2(a, b, n);  % Cheby2 mesh

    [ky_uni, k1_uni, Lam_uni, H_uni, ~] = run_conditioning(x_uni, x_eval);  % Uniform
    [ky_ch1, k1_ch1, Lam_ch1, H_ch1, ~] = run_conditioning(x_ch1, x_eval);  % Cheby1
    [ky_ch2, k1_ch2, Lam_ch2, H_ch2, ~] = run_conditioning(x_ch2, x_eval);  % Cheby2

    Lam_uni_all(ni) = Lam_uni;  H_uni_all(ni) = H_uni;  % Store uniform stats
    Lam_ch1_all(ni) = Lam_ch1;  H_ch1_all(ni) = H_ch1;  % Store Cheby1 stats
    Lam_ch2_all(ni) = Lam_ch2;  H_ch2_all(ni) = H_ch2;  % Store Cheby2 stats

    fprintf('%-6d %-10s %-14.4e %-14.4e %-14.4e %-14.4e\n', n, 'Uniform', ...  % Uniform row
        Lam_uni, H_uni, max(ky_uni(isfinite(ky_uni))), max(k1_uni));
    fprintf('%-6d %-10s %-14.4e %-14.4e %-14.4e %-14.4e\n', n, 'Cheby1',  ...  % Cheby1 row
        Lam_ch1, H_ch1, max(ky_ch1(isfinite(ky_ch1))), max(k1_ch1));
    fprintf('%-6d %-10s %-14.4e %-14.4e %-14.4e %-14.4e\n', n, 'Cheby2',  ...  % Cheby2 row
        Lam_ch2, H_ch2, max(ky_ch2(isfinite(ky_ch2))), max(k1_ch2));
    fprintf('%s\n', repmat('-', 1, 76));  % Separator between degrees
end
fprintf('\n');  % Blank line

% Figure 1: Lambda_n vs n for all meshes
figure(1); clf; hold on; grid on;                                                                       % Create figure
semilogy(n_sweep, Lam_uni_all, '-o', 'LineWidth', 1.5, 'DisplayName', 'Uniform');                      % Uniform
semilogy(n_sweep, Lam_ch1_all, '-s', 'LineWidth', 1.5, 'DisplayName', 'Cheby1');                       % Cheby1
semilogy(n_sweep, Lam_ch2_all, '-d', 'LineWidth', 1.5, 'DisplayName', 'Cheby2');                       % Cheby2
xticks(n_sweep);                                                                                        % Set x ticks
xlabel('n'); ylabel('\Lambda_n');                                                                       % Label axes
title('Figure 1: Lebesgue constant \Lambda_n vs n, f4');                                               % Set title
legend('Location','best'); hold off;                                                                    % Add legend

% Figure 2: H_n vs n for all meshes
figure(2); clf; hold on; grid on;                                                                       % Create figure
semilogy(n_sweep, H_uni_all, '-o', 'LineWidth', 1.5, 'DisplayName', 'Uniform');                        % Uniform
semilogy(n_sweep, H_ch1_all, '-s', 'LineWidth', 1.5, 'DisplayName', 'Cheby1');                         % Cheby1
semilogy(n_sweep, H_ch2_all, '-d', 'LineWidth', 1.5, 'DisplayName', 'Cheby2');                         % Cheby2
xticks(n_sweep);                                                                                        % Set x ticks
xlabel('n'); ylabel('H_n');                                                                             % Label axes
title('Figure 2: H_n norm vs n, f4');                                                                  % Set title
legend('Location','best'); hold off;                                                                    % Add legend

% ============================================================
% SUBTASK 2: Accuracy and stability at n=29 (Higham-style plots)
% ============================================================

fprintf('------------------------------------------------------------\n');   % Print separator
fprintf('SUBTASK 2: Accuracy and stability at n=%d\n', n_higham);            % Section header
fprintf('------------------------------------------------------------\n\n'); % Print separator

x_uni_h = mesh_uniform(a, b, n_higham);  % Uniform mesh at n=29
x_ch1_h = mesh_cheby1(a, b, n_higham);  % Cheby1 mesh at n=29
x_ch2_h = mesh_cheby2(a, b, n_higham);  % Cheby2 mesh at n=29

% Stability bound for Bary Form 2 (Set 8 Section 5.4 / Higham 2004):
%   (3n+4)*u*kappa_y + (3n+2)*u*kappa_1, multiplied by |f_exact|
u = double(u_sp);
[ky_uni_h, k1_uni_h, ~, ~, ~] = run_conditioning(x_uni_h, x_eval);  % Uniform conditioning at n=29
[ky_ch1_h, k1_ch1_h, ~, ~, ~] = run_conditioning(x_ch1_h, x_eval);  % Cheby1 conditioning at n=29

stab_uni_h = ((3*n_higham+4)*u*ky_uni_h + (3*n_higham+2)*u*k1_uni_h) .* abs(f_exact);  % Bary2 bound uniform
stab_ch1_h = ((3*n_higham+4)*u*ky_ch1_h + (3*n_higham+2)*u*k1_ch1_h) .* abs(f_exact);  % Bary2 bound Cheby1

[p_b2_uni_h, p_inc_uni_h, p_dec_uni_h, p_lej_uni_h] = eval_all_forms(x_uni_h, x_eval, 'uniform');  % Uniform forms
[p_b2_ch1_h, p_inc_ch1_h, p_dec_ch1_h, p_lej_ch1_h] = eval_all_forms(x_ch1_h, x_eval, 'cheby1');  % Cheby1 forms
[p_b2_ch2_h, p_inc_ch2_h, p_dec_ch2_h, p_lej_ch2_h] = eval_all_forms(x_ch2_h, x_eval, 'cheby2');  % Cheby2 forms

fprintf('Table 2: max |p(x) - f4(x)| for all forms at n=%d\n', n_higham);                        % Table header
fprintf('%-24s %-14s %-14s %-14s\n', 'Form', 'Uniform', 'Cheby1', 'Cheby2');                     % Column headers
fprintf('%s\n', repmat('-', 1, 68));                                                               % Separator
fprintf('%-24s %-14.4e %-14.4e %-14.4e\n', 'Bary Form 2',       ...                              % Bary2 row
    max(abs(double(p_b2_uni_h)  - f_exact)), ...
    max(abs(double(p_b2_ch1_h)  - f_exact)), ...
    max(abs(double(p_b2_ch2_h)  - f_exact)));
fprintf('%-24s %-14.4e %-14.4e %-14.4e\n', 'Newton increasing', ...                              % Newton inc row
    max(abs(double(p_inc_uni_h) - f_exact)), ...
    max(abs(double(p_inc_ch1_h) - f_exact)), ...
    max(abs(double(p_inc_ch2_h) - f_exact)));
fprintf('%-24s %-14.4e %-14.4e %-14.4e\n', 'Newton decreasing', ...                              % Newton dec row
    max(abs(double(p_dec_uni_h) - f_exact)), ...
    max(abs(double(p_dec_ch1_h) - f_exact)), ...
    max(abs(double(p_dec_ch2_h) - f_exact)));
fprintf('%-24s %-14.4e %-14.4e %-14.4e\n', 'Newton Leja',       ...                              % Newton Leja row
    max(abs(double(p_lej_uni_h) - f_exact)), ...
    max(abs(double(p_lej_ch1_h) - f_exact)), ...
    max(abs(double(p_lej_ch2_h) - f_exact)));
fprintf('\n');  % Blank line

% Figure 3: Higham-style uniform n=29
figure(3); clf; hold on; grid on;                                                                  % Create figure
semilogy(x_eval, max(abs(double(p_b2_uni_h)  - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Bary2');       % Bary2
semilogy(x_eval, max(abs(double(p_inc_uni_h) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton inc');  % Newton inc
semilogy(x_eval, max(abs(double(p_dec_uni_h) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton dec');  % Newton dec
semilogy(x_eval, max(abs(double(p_lej_uni_h) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton Leja'); % Newton Leja
semilogy(x_eval, max(stab_uni_h, 1e-20), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Stability bound');             % Stability bound
xlabel('x'); ylabel('|p(x) - f4(x)|');                                                            % Label axes
title(sprintf('Figure 3: Pointwise error, f4, uniform mesh, n=%d', n_higham));                    % Set title
legend('Location','best'); hold off;                                                               % Add legend

% Figure 4: Higham-style Cheby1 n=29
figure(4); clf; hold on; grid on;                                                                  % Create figure
semilogy(x_eval, max(abs(double(p_b2_ch1_h)  - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Bary2');       % Bary2
semilogy(x_eval, max(abs(double(p_inc_ch1_h) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton inc');  % Newton inc
semilogy(x_eval, max(abs(double(p_dec_ch1_h) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton dec');  % Newton dec
semilogy(x_eval, max(abs(double(p_lej_ch1_h) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton Leja'); % Newton Leja
semilogy(x_eval, max(stab_ch1_h, 1e-20), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Stability bound');             % Stability bound
xlabel('x'); ylabel('|p(x) - f4(x)|');                                                            % Label axes
title(sprintf('Figure 4: Pointwise error, f4, Cheby1 mesh, n=%d', n_higham));                     % Set title
legend('Location','best'); hold off;                                                               % Add legend

% ============================================================
% SUBTASK 3: Convergence investigation
% ============================================================

fprintf('------------------------------------------------------------\n');   % Print separator
fprintf('SUBTASK 3: Convergence investigation\n');                            % Section header
fprintf('------------------------------------------------------------\n\n'); % Print separator

% Use a finer n sweep for convergence investigation
n_conv = [5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 75, 100];  % n values for convergence sweep

inf_err_uni  = zeros(length(n_conv), 1);  % Inf-norm error storage uniform
inf_err_ch1  = zeros(length(n_conv), 1);  % Inf-norm error storage Cheby1
inf_err_ch2  = zeros(length(n_conv), 1);  % Inf-norm error storage Cheby2

for ni = 1 : length(n_conv)  % Loop over each degree
    n = n_conv(ni);           % Get current degree

    x_uni_c = mesh_uniform(a, b, n);  % Uniform mesh
    x_ch1_c = mesh_cheby1(a, b, n);  % Cheby1 mesh
    x_ch2_c = mesh_cheby2(a, b, n);  % Cheby2 mesh

    x_s_uni = single(x_uni_c);  fv_uni_c = f4(x_s_uni);  % Uniform nodes and values
    x_s_ch1 = single(x_ch1_c);  fv_ch1_c = f4(x_s_ch1);  % Cheby1 nodes and values
    x_s_ch2 = single(x_ch2_c);  fv_ch2_c = f4(x_s_ch2);  % Cheby2 nodes and values

    b_uni = bary2_coeffs(double(x_uni_c), n, 'uniform');  % Bary2 weights uniform O(n)
    b_ch1 = bary2_coeffs(double(x_ch1_c), n, 'cheby1');   % Bary2 weights Cheby1 O(n)
    b_ch2 = bary2_coeffs(double(x_ch2_c), n, 'cheby2');   % Bary2 weights Cheby2 O(n)

    p_uni_c = bary2_eval(x_s_uni, fv_uni_c, single(b_uni), single(x_eval));  % Bary2 uniform
    p_ch1_c = bary2_eval(x_s_ch1, fv_ch1_c, single(b_ch1), single(x_eval));  % Bary2 Cheby1
    p_ch2_c = bary2_eval(x_s_ch2, fv_ch2_c, single(b_ch2), single(x_eval));  % Bary2 Cheby2

    inf_err_uni(ni) = max(abs(double(p_uni_c) - f_exact));  % Inf-norm error uniform
    inf_err_ch1(ni) = max(abs(double(p_ch1_c) - f_exact));  % Inf-norm error Cheby1
    inf_err_ch2(ni) = max(abs(double(p_ch2_c) - f_exact));  % Inf-norm error Cheby2
end

% Print convergence table for Cheby1 and Cheby2 only (uniform diverges)
fprintf('Table 3: ||f4 - p_n||_inf vs n (Bary2, convergent meshes)\n');              % Table header
fprintf('%-6s %-18s %-18s %-18s\n', 'n', 'Uniform', 'Cheby1', 'Cheby2');             % Column headers
fprintf('%s\n', repmat('-', 1, 62));                                                  % Separator
for ni = 1 : length(n_conv)  % Loop over each degree
    fprintf('%-6d %-18.4e %-18.4e %-18.4e\n', n_conv(ni), ...                        % Print row
        inf_err_uni(ni), inf_err_ch1(ni), inf_err_ch2(ni));
end
fprintf('\n');  % Blank line

% Report n needed to reach each target accuracy level for Cheby1
fprintf('Cheby1 - n required to reach target accuracy levels:\n');  % Print header
for ti = 1 : length(target_levels)  % Loop over each target level
    idx = find(inf_err_ch1 <= target_levels(ti), 1, 'first');  % Find first n meeting target
    if ~isempty(idx)  % If a qualifying n was found
        fprintf('  ||f4-p_n||_inf <= %.0e  ->  n = %d  (error = %.4e)\n', ...  % Print result
            target_levels(ti), n_conv(idx), inf_err_ch1(idx));
    else  % Target not reached within sweep range
        fprintf('  ||f4-p_n||_inf <= %.0e  ->  not reached in sweep\n', target_levels(ti));  % Print result
    end
end
fprintf('\n');  % Blank line

% Single precision floor: the Bary2 error cannot go below ~ u * kappa(x,n,y) * |f(x)|
% (Set 8 Section 5.4). For well-conditioned Chebyshev meshes kappa ~ 1, so the
% floor approaches u ~ 6e-8, but will be higher wherever kappa_y is large.
floor_ch1 = min(inf_err_ch1);  % Minimum observed error Cheby1
floor_ch2 = min(inf_err_ch2);  % Minimum observed error Cheby2
fprintf('Single precision floor observations (threshold ~ u * kappa(x,n,y)):\n');  % Print header
fprintf('  Min ||f4-p_n||_inf Cheby1 = %.4e  (at n=%d)\n', ...                     % Cheby1 floor
    floor_ch1, n_conv(inf_err_ch1 == floor_ch1));
fprintf('  Min ||f4-p_n||_inf Cheby2 = %.4e  (at n=%d)\n\n', ...                   % Cheby2 floor
    floor_ch2, n_conv(inf_err_ch2 == floor_ch2));
fprintf('  Single precision u = %.4e\n\n', double(u_sp));  % Print unit roundoff

% Figure 5: inf-norm error vs n all meshes
figure(5); clf; hold on; grid on;                                                                        % Create figure
semilogy(n_conv, inf_err_uni, '-o', 'LineWidth', 1.5, 'DisplayName', 'Uniform (Runge divergence)');     % Uniform
semilogy(n_conv, inf_err_ch1, '-s', 'LineWidth', 1.5, 'DisplayName', 'Cheby1 (convergent)');            % Cheby1
semilogy(n_conv, inf_err_ch2, '-d', 'LineWidth', 1.5, 'DisplayName', 'Cheby2 (convergent)');            % Cheby2
yline(double(u_sp), 'k--', 'LineWidth', 1.2, 'DisplayName', 'Single precision u');                      % Precision floor reference
xlabel('n'); ylabel('||f4 - p_n||_\infty');                                                              % Label axes
title('Figure 5: Convergence of Bary2 to f4 as n increases');                                           % Set title
legend('Location','best'); hold off;                                                                     % Add legend

% Figure 6: Convergence zoom - Cheby1 and Cheby2 only with target levels
figure(6); clf; hold on; grid on;                                                                        % Create figure
semilogy(n_conv, inf_err_ch1, '-s', 'LineWidth', 1.5, 'DisplayName', 'Cheby1');                         % Cheby1
semilogy(n_conv, inf_err_ch2, '-d', 'LineWidth', 1.5, 'DisplayName', 'Cheby2');                         % Cheby2
yline(double(u_sp), 'k--', 'LineWidth', 1.2, 'DisplayName', 'Single precision u (floor)');              % Precision floor
for ti = 1 : length(target_levels)  % Loop over target levels
    yline(target_levels(ti), ':', 'LineWidth', 1.0, ...                                                  % Draw target line
        'DisplayName', sprintf('Target %.0e', target_levels(ti)));
end
xlabel('n'); ylabel('||f4 - p_n||_\infty');                                                              % Label axes
title('Figure 6: Convergence zoom - Cheby meshes and precision floor');                                  % Set title
legend('Location','best'); hold off;                                                                     % Add legend

fprintf('DONE. Produced Tables 1-3 and Figures 1-6.\n');  % Print completion message
