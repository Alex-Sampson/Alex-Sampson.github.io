% task1driver_p2.m
%
% Task 1: Code design, complexity, and empirical validation.
% Validates all interpolation subroutines on known polynomials before
% the main experiments in Tasks 2-5.
%
% Subroutines exercised:
%   mesh_uniform(a, b, n)
%   mesh_cheby1(a, b, n)
%   mesh_cheby2(a, b, n)
%   leja_order(x)
%   product_eval(x, alpha, rho)
%   bary2_eval(x_nodes, f_vals, beta, x)
%   newton_coeffs(x_nodes, f_vals)
%   newton_eval(x_nodes, c, x)
%   bary1_eval_dp(x_nodes, f_vals, gamma, x)
%   cond_kappa(x_nodes, f_vals, gamma, x, denom)
%   lebesgue_lambda(kappa_1)
%   hilbert_H(kappa_y)
%
% Output:
%   Table 1: Mesh point validation (uniform, Cheby1, Cheby2, Leja)
%   Table 2: product_eval vs exact for f1 and f2
%   Table 3: bary2_eval, newton_eval accuracy vs exact for f1
%   Table 4: bary1_eval_dp accuracy vs exact for f1
%   Figure 1: Mesh point distributions for n=10
%   Figure 2: Error |p(x) - f(x)| for all forms on f1, uniform mesh
%   Figure 3: Error |p(x) - f(x)| for all forms on f1, Cheby1 mesh
%   Figure 4: Lebesgue function kappa(x,n,1) for uniform vs Cheby1

clear; clc; close all;  % Clear workspace, command window, and close all figures

% ----------------  settings  ----------------
a     = 0.5;   % Left endpoint of interpolation interval (contains rho=2 for f1)
b     = 3.5;   % Right endpoint of interpolation interval
n     = 9;     % Degree of interpolation for f1 (10 points)
rho   = 2.0;   % Root parameter for f1
d     = 9;     % Degree parameter for f1
n_pts = 500;   % Number of query points for evaluation and plotting

x_eval = linspace(a, b, n_pts)';  % Dense query grid over [a, b]

fprintf('\n============================================================\n');  % Print separator
fprintf('TASK 1: Code Validation\n');                                          % Print task header
fprintf('a=%.4g, b=%.4g, n=%d, rho=%.4g, d=%d\n', a, b, n, rho, d);          % Print settings
fprintf('============================================================\n\n');   % Print separator

% ============================================================
% VALIDATION 1: Mesh generation
% ============================================================

fprintf('------------------------------------------------------------\n');  % Print separator
fprintf('VALIDATION 1: Mesh generation\n');                                  % Print section header
fprintf('------------------------------------------------------------\n\n'); % Print separator

x_uni   = mesh_uniform(a, b, n);   % Generate uniform mesh
x_ch1   = mesh_cheby1(a, b, n);    % Generate Chebyshev first kind mesh
x_ch2   = mesh_cheby2(a, b, n);    % Generate Chebyshev second kind mesh
x_leja  = leja_order(x_uni);       % Reorder uniform mesh into Leja ordering

fprintf('Table 1: Mesh points for n=%d on [%.4g, %.4g]\n', n, a, b);  % Print table header
fprintf('%-6s %-14s %-14s %-14s %-14s\n', 'i', 'Uniform', 'Cheby1', 'Cheby2', 'Leja');  % Print column headers
fprintf('%s\n', repmat('-',1,62));  % Print separator line
for i = 1 : n+1  % Loop over each mesh point
    fprintf('%-6d %-14.8f %-14.8f %-14.8f %-14.8f\n', ...  % Print row of mesh values
        i, x_uni(i), x_ch1(i), x_ch2(i), x_leja(i));
end
fprintf('\n');  % Blank line after table

% Check endpoints are correct for all meshes
assert(abs(x_uni(1)  - a) < 1e-12, 'Uniform: left endpoint wrong');   % Validate uniform left endpoint
assert(abs(x_uni(end)- b) < 1e-12, 'Uniform: right endpoint wrong');  % Validate uniform right endpoint
assert(abs(x_ch2(1)  - a) < 1e-12, 'Cheby2: left endpoint wrong');    % Cheby2 must include endpoints
assert(abs(x_ch2(end)- b) < 1e-12, 'Cheby2: right endpoint wrong');   % Cheby2 must include endpoints
fprintf('Endpoint checks passed for uniform and Cheby2 meshes.\n\n');  % Confirm checks passed

% Figure 1: mesh point distributions
figure(1); clf; hold on; grid on;  % Create figure and enable grid
plot(x_uni,  ones(n+1,1)*1, 'o', 'LineWidth', 1.5, 'DisplayName', 'Uniform');   % Plot uniform points
plot(x_ch1,  ones(n+1,1)*2, 's', 'LineWidth', 1.5, 'DisplayName', 'Cheby1');    % Plot Cheby1 points
plot(x_ch2,  ones(n+1,1)*3, 'd', 'LineWidth', 1.5, 'DisplayName', 'Cheby2');    % Plot Cheby2 points
plot(x_leja, ones(n+1,1)*4, '^', 'LineWidth', 1.5, 'DisplayName', 'Leja');      % Plot Leja points
legend('Location','best');  % Add legend
xlabel('x'); yticks([1 2 3 4]); yticklabels({'Uniform','Cheby1','Cheby2','Leja'});  % Label axes
title('Figure 1: Mesh point distributions (n=9)');  % Set title
hold off;  % Release hold

% ============================================================
% VALIDATION 2: product_eval vs exact for f1 and f2
% ============================================================

fprintf('------------------------------------------------------------\n');  % Print separator
fprintf('VALIDATION 2: product_eval accuracy\n');                            % Print section header
fprintf('------------------------------------------------------------\n\n'); % Print separator

% f1: (x - rho)^d, exact via MATLAB built-in in double
p_prod_f1  = product_eval(double(x_eval), 1.0, rho*ones(d,1));  % Evaluate f1 via product form
p_exact_f1 = (double(x_eval) - rho) .^ d;                       % Evaluate f1 exactly in double

% f2: prod_{i=1}^{d}(x-i), exact via product_eval itself in double (known accurate)
p_prod_f2  = product_eval(double(x_eval), 1.0, (1:d)');  % Evaluate f2 via product form
p_exact_f2 = ones(n_pts, 1);                              % Initialize exact f2 values
for i = 1 : d                                             % Loop over each root
    p_exact_f2 = p_exact_f2 .* (double(x_eval) - i);     % Accumulate product in double
end

err_f1 = max(abs(p_prod_f1 - p_exact_f1));  % Max absolute error for f1
err_f2 = max(abs(p_prod_f2 - p_exact_f2));  % Max absolute error for f2

fprintf('Table 2: product_eval max error over %d query points\n', n_pts);  % Print table header
fprintf('%-10s %-20s\n', 'Function', 'max |product - exact|');              % Print column headers
fprintf('%s\n', repmat('-',1,32));                                           % Print separator
fprintf('%-10s %-20.6e\n', 'f1', err_f1);  % Print f1 error
fprintf('%-10s %-20.6e\n', 'f2', err_f2);  % Print f2 error
fprintf('\n');  % Blank line

% ============================================================
% VALIDATION 3: Interpolation accuracy - all forms on f1
% ============================================================

fprintf('------------------------------------------------------------\n');  % Print separator
fprintf('VALIDATION 3: Interpolation accuracy on f1 (all forms)\n');        % Print section header
fprintf('------------------------------------------------------------\n\n'); % Print separator

% Reference exact values (double precision product form)
f_exact_dp = product_eval(double(x_eval), 1.0, rho*ones(d,1));  % Exact f1 values at query points

% --- Uniform mesh ---
x_uni_s  = single(x_uni);                                % Cast uniform nodes to single
fv_uni   = single(f1(x_uni_s, rho, d));                  % Function values at nodes (single)
gamma_uni = bary1_coeffs(double(x_uni), n, @(x) (x-rho).^d); % Compute Bary1 weights (double)
beta_uni  = bary2_coeffs(double(x_uni), n, 'uniform');         % Compute Bary2 weights O(n)

% Barycentric Form 2 on uniform mesh
p_b2_uni = bary2_eval(x_uni_s, fv_uni, single(beta_uni), single(x_eval));  % Evaluate Bary2

% Newton increasing order on uniform mesh
c_inc    = newton_coeffs(x_uni_s, fv_uni);                          % Compute Newton coefficients
p_n_inc  = newton_eval(x_uni_s, c_inc, single(x_eval));             % Evaluate Newton increasing

% Newton decreasing order on uniform mesh
x_dec    = flipud(x_uni_s);                                          % Reverse node order
fv_dec   = flipud(fv_uni);                                           % Reverse function values
c_dec    = newton_coeffs(x_dec, fv_dec);                            % Compute Newton coefficients
p_n_dec  = newton_eval(x_dec, c_dec, single(x_eval));               % Evaluate Newton decreasing

% Newton Leja order on uniform mesh
x_lej_s  = single(x_leja);                                          % Cast Leja nodes to single
fv_lej   = single(f1(x_lej_s, rho, d));                             % Function values at Leja nodes
c_lej    = newton_coeffs(x_lej_s, fv_lej);                          % Compute Newton coefficients
p_n_lej  = newton_eval(x_lej_s, c_lej, single(x_eval));             % Evaluate Newton Leja

err_b2_uni  = max(abs(double(p_b2_uni) - f_exact_dp));   % Max error Bary2 uniform
err_inc_uni = max(abs(double(p_n_inc)  - f_exact_dp));   % Max error Newton increasing uniform
err_dec_uni = max(abs(double(p_n_dec)  - f_exact_dp));   % Max error Newton decreasing uniform
err_lej_uni = max(abs(double(p_n_lej)  - f_exact_dp));   % Max error Newton Leja uniform

fprintf('Table 3: max |p(x) - f1(x)| over %d query points (single precision forms)\n', n_pts);  % Table header
fprintf('%-28s %-14s %-14s\n', 'Form', 'Uniform mesh', 'Cheby1 mesh');  % Column headers
fprintf('%s\n', repmat('-',1,58));  % Separator

% --- Cheby1 mesh ---
x_ch1_s  = single(x_ch1);                                % Cast Cheby1 nodes to single
fv_ch1   = single(f1(x_ch1_s, rho, d));                  % Function values at Cheby1 nodes
gamma_ch1 = bary1_coeffs(double(x_ch1), n, @(x) (x-rho).^d); % Compute Bary1 weights (double)
beta_ch1  = bary2_coeffs(double(x_ch1), n, 'cheby1');          % Compute Bary2 weights O(n)

p_b2_ch1 = bary2_eval(x_ch1_s, fv_ch1, single(beta_ch1), single(x_eval));  % Bary2 Cheby1

x_dec_ch1  = flipud(x_ch1_s);                            % Reverse Cheby1 node order
fv_dec_ch1 = flipud(fv_ch1);                             % Reverse function values
c_inc_ch1  = newton_coeffs(x_ch1_s, fv_ch1);            % Newton coefficients increasing
c_dec_ch1  = newton_coeffs(x_dec_ch1, fv_dec_ch1);      % Newton coefficients decreasing
x_lej_ch1  = single(leja_order(x_ch1));                  % Leja order of Cheby1 nodes
fv_lej_ch1 = single(f1(x_lej_ch1, rho, d));             % Function values at Leja-ordered nodes
c_lej_ch1  = newton_coeffs(x_lej_ch1, fv_lej_ch1);      % Newton coefficients Leja

p_n_inc_ch1 = newton_eval(x_ch1_s,    c_inc_ch1,  single(x_eval));  % Newton increasing Cheby1
p_n_dec_ch1 = newton_eval(x_dec_ch1,  c_dec_ch1,  single(x_eval));  % Newton decreasing Cheby1
p_n_lej_ch1 = newton_eval(x_lej_ch1,  c_lej_ch1,  single(x_eval));  % Newton Leja Cheby1

err_b2_ch1  = max(abs(double(p_b2_ch1)     - f_exact_dp));  % Max error Bary2 Cheby1
err_inc_ch1 = max(abs(double(p_n_inc_ch1)  - f_exact_dp));  % Max error Newton increasing Cheby1
err_dec_ch1 = max(abs(double(p_n_dec_ch1)  - f_exact_dp));  % Max error Newton decreasing Cheby1
err_lej_ch1 = max(abs(double(p_n_lej_ch1)  - f_exact_dp));  % Max error Newton Leja Cheby1

fprintf('%-28s %-14.4e %-14.4e\n', 'Barycentric Form 2',    err_b2_uni,  err_b2_ch1);   % Print Bary2 errors
fprintf('%-28s %-14.4e %-14.4e\n', 'Newton (increasing)',   err_inc_uni, err_inc_ch1);  % Print Newton inc errors
fprintf('%-28s %-14.4e %-14.4e\n', 'Newton (decreasing)',   err_dec_uni, err_dec_ch1);  % Print Newton dec errors
fprintf('%-28s %-14.4e %-14.4e\n', 'Newton (Leja)',         err_lej_uni, err_lej_ch1);  % Print Newton Leja errors
fprintf('\n');  % Blank line

% ============================================================
% VALIDATION 4: bary1_eval_dp accuracy and kappa outputs
% ============================================================

fprintf('------------------------------------------------------------\n');  % Print separator
fprintf('VALIDATION 4: bary1_eval_dp and condition number outputs\n');       % Print section header
fprintf('------------------------------------------------------------\n\n'); % Print separator

fv_uni_dp  = (double(x_uni) - rho) .^ d;               % Double precision function values at nodes
[p_dp, numer_abs, denom_b1] = bary1_eval_dp( ...        % Evaluate Bary1 in double
    double(x_uni), fv_uni_dp, gamma_uni, double(x_eval));

denom_prod = product_eval(double(x_eval), 1.0, rho*ones(d,1));  % Product form denominator

err_dp = max(abs(p_dp - f_exact_dp));  % Max error of bary1_eval_dp vs exact

fprintf('Table 4: bary1_eval_dp max error vs product_eval reference\n');  % Print table header
fprintf('  max |bary1_dp - exact| = %.6e\n\n', err_dp);                   % Print error

% Condition numbers using product form denominator
[kappa_y, kappa_1] = cond_kappa( ...                                  % Compute condition numbers
    double(x_uni), fv_uni_dp, gamma_uni, double(x_eval), abs(denom_prod));

Lambda_n = lebesgue_lambda(kappa_1);  % Compute Lebesgue constant
H_n      = hilbert_H(kappa_y);        % Compute H_n norm

fprintf('  Lebesgue constant Lambda_n (uniform, n=%d) = %.6e\n', n, Lambda_n);  % Print Lambda_n
fprintf('  H_n norm          H_n      (uniform, n=%d) = %.6e\n\n', n, H_n);     % Print H_n

% ============================================================
% FIGURES 2-4: Error and Lebesgue function plots
% ============================================================

% Figure 2: pointwise error on uniform mesh
figure(2); clf; hold on; grid on;  % Create figure
semilogy(x_eval, max(abs(double(p_b2_uni) - f_exact_dp), 1e-20),  'LineWidth', 1.5, 'DisplayName', 'Bary2');         % Plot Bary2 error
semilogy(x_eval, max(abs(double(p_n_inc)  - f_exact_dp), 1e-20),  'LineWidth', 1.5, 'DisplayName', 'Newton inc');    % Plot Newton inc error
semilogy(x_eval, max(abs(double(p_n_dec)  - f_exact_dp), 1e-20),  'LineWidth', 1.5, 'DisplayName', 'Newton dec');    % Plot Newton dec error
semilogy(x_eval, max(abs(double(p_n_lej)  - f_exact_dp), 1e-20),  'LineWidth', 1.5, 'DisplayName', 'Newton Leja');   % Plot Newton Leja error
xlabel('x'); ylabel('|p(x) - f1(x)|');                 % Label axes
title('Figure 2: Pointwise error, f1, uniform mesh');  % Set title
legend('Location','best'); hold off;                    % Add legend

% Figure 3: pointwise error on Cheby1 mesh
figure(3); clf; hold on; grid on;  % Create figure
semilogy(x_eval, max(abs(double(p_b2_ch1)    - f_exact_dp), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Bary2');        % Plot Bary2 error
semilogy(x_eval, max(abs(double(p_n_inc_ch1) - f_exact_dp), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton inc');   % Plot Newton inc error
semilogy(x_eval, max(abs(double(p_n_dec_ch1) - f_exact_dp), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton dec');   % Plot Newton dec error
semilogy(x_eval, max(abs(double(p_n_lej_ch1) - f_exact_dp), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton Leja');  % Plot Newton Leja error
xlabel('x'); ylabel('|p(x) - f1(x)|');                  % Label axes
title('Figure 3: Pointwise error, f1, Cheby1 mesh');    % Set title
legend('Location','best'); hold off;                     % Add legend

% Figure 4: Lebesgue function for uniform vs Cheby1
[~, kappa_1_ch1] = cond_kappa( ...                                  % Condition numbers for Cheby1
    double(x_ch1), (double(x_ch1)-rho).^d, gamma_ch1, ...
    double(x_eval), abs(product_eval(double(x_eval), 1.0, rho*ones(d,1))));

figure(4); clf; hold on; grid on;  % Create figure
plot(x_eval, kappa_1,     'LineWidth', 1.5, 'DisplayName', 'Uniform');   % Plot Lebesgue function uniform
plot(x_eval, kappa_1_ch1, 'LineWidth', 1.5, 'DisplayName', 'Cheby1');    % Plot Lebesgue function Cheby1
xlabel('x'); ylabel('\kappa(x,n,1)');                        % Label axes
title('Figure 4: Lebesgue function, uniform vs Cheby1');     % Set title
legend('Location','best'); hold off;                          % Add legend

fprintf('DONE. Produced Tables 1-4 and Figures 1-4.\n');  % Print completion message
