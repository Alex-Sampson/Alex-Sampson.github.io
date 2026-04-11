% task3driver_p2.m
%
% Task 3: Interpolation of f2(x; d) = prod_{i=1}^{d} (x - i)
%
% Subroutines used:
%   mesh_uniform, mesh_cheby1, mesh_cheby2, leja_order
%   product_eval, bary1_coeffs, bary1_eval_dp
%   bary2_eval, newton_coeffs, newton_eval
%   cond_kappa, lebesgue_lambda, hilbert_H
%   f2
%
% Subtask 1: Interpolation on uniform, Cheby1, Cheby2 meshes for multiple
%            degrees (d=10, 21, 25) including overdetermined tests
% Subtask 2: Conditioning - kappa(x,n,y), kappa(x,n,1), Lambda_n, H_n
% Subtask 3: Accuracy and stability of single precision forms vs exact
%
% Output:
%   Table 1: Lambda_n and H_n for all degrees and mesh types
%   Table 2: Max error for all forms, all degrees, all mesh types
%   Table 3: Overdetermined test errors (more points than degree+1)
%   Figure 1: Lebesgue function, uniform vs Cheby1, d=25
%   Figure 2: kappa(x,n,y), uniform vs Cheby1, d=25
%   Figure 3: Pointwise error all forms, uniform mesh, d=25 (Higham-style)
%   Figure 4: Pointwise error all forms, Cheby1 mesh, d=25 (Higham-style)
%   Figure 5: Lambda_n vs degree for uniform and Cheby1
%   Figure 6: Max error vs degree for Bary2 on uniform and Cheby1

clear; clc; close all;  % Clear workspace, command window, and close all figures

% ----------------  settings  ----------------
d_list      = [10, 21, 25];   % Degrees to test (at least two > 20 as required)
d_over_list = [15, 30];       % Overdetermined: more points than degree+1
d_over_deg  = 10;             % True degree of f2 for overdetermined tests
a           = 0.0;            % Left endpoint (contains roots 1..d)
b_func      = @(d) d + 1.0;  % Right endpoint as function of d (contains all roots)
n_pts       = 1000;           % Number of query points
u_sp        = single(eps('single'));  % Unit roundoff for single precision

fprintf('\n============================================================\n');  % Print separator
fprintf('TASK 3: f2(x; d) = prod_{i=1}^{d} (x - i)\n');                      % Print task header
fprintf('Degrees tested: ');  fprintf('%d ', d_list);  fprintf('\n');          % Print degree list
fprintf('============================================================\n\n');   % Print separator

% ============================================================
% LOCAL HELPERS
% ============================================================

function [p_b2, p_inc, p_dec, p_lej] = eval_all_forms(x_nodes, d, x_eval, mesh_type)
    % Evaluate all four single-precision forms for f2 on a given mesh.
    x_s    = single(x_nodes);                              % Cast nodes to single
    fv     = f2(x_s, d);                                   % Function values at nodes (single)
    beta   = bary2_coeffs(double(x_nodes), length(x_nodes)-1, mesh_type); % Bary2 weights O(n)

    p_b2   = bary2_eval(x_s, fv, single(beta), single(x_eval));  % Barycentric Form 2

    c_inc  = newton_coeffs(x_s, fv);                               % Newton coeffs increasing
    p_inc  = newton_eval(x_s, c_inc, single(x_eval));              % Newton increasing

    x_dec  = flipud(x_s);  fv_dec = flipud(fv);                    % Reverse order
    c_dec  = newton_coeffs(x_dec, fv_dec);                         % Newton coeffs decreasing
    p_dec  = newton_eval(x_dec, c_dec, single(x_eval));            % Newton decreasing

    x_lej  = single(leja_order(x_nodes));                          % Leja ordering
    fv_lej = f2(x_lej, d);                                         % Function values at Leja nodes
    c_lej  = newton_coeffs(x_lej, fv_lej);                         % Newton coeffs Leja
    p_lej  = newton_eval(x_lej, c_lej, single(x_eval));            % Newton Leja
end

function [kappa_y, kappa_1, Lambda_n, H_n] = run_conditioning(x_nodes, d, x_eval)
    % Compute condition numbers and summary stats for f2 on a given mesh.
    fv_dp  = double(f2(double(x_nodes), d));                        % Function values in double
    gamma  = bary1_coeffs(double(x_nodes), length(x_nodes)-1, ...
                          @(t) f2(t, d));                           % Bary1 weights (double)
    denom  = abs(product_eval(double(x_eval), 1.0, (1:d)'));        % Exact product form denominator

    [kappa_y, kappa_1] = cond_kappa(double(x_nodes), fv_dp, gamma, ...
                                     double(x_eval), denom);        % Condition numbers
    Lambda_n = lebesgue_lambda(kappa_1);  % Lebesgue constant
    H_n      = hilbert_H(kappa_y);        % H_n norm
end

% ============================================================
% SUBTASKS 1 & 2: Multi-degree conditioning sweep
% ============================================================

fprintf('------------------------------------------------------------\n');   % Print separator
fprintf('SUBTASKS 1 & 2: Conditioning sweep over degrees\n');                % Section header
fprintf('------------------------------------------------------------\n\n'); % Print separator

fprintf('Table 1: Conditioning summary for f2 (all degrees and meshes)\n');                      % Table header
fprintf('%-6s %-10s %-14s %-14s %-14s %-14s\n', ...                                              % Column headers
    'd', 'Mesh', 'Lambda_n', 'H_n', 'max kappa_y', 'max kappa_1');
fprintf('%s\n', repmat('-', 1, 76));  % Separator

Lam_uni_all = zeros(length(d_list), 1);  % Store Lambda_n uniform for each degree
Lam_ch1_all = zeros(length(d_list), 1);  % Store Lambda_n Cheby1 for each degree
err_b2_uni_all  = zeros(length(d_list), 1);  % Store Bary2 uniform errors
err_b2_ch1_all  = zeros(length(d_list), 1);  % Store Bary2 Cheby1 errors

for di = 1 : length(d_list)  % Loop over each degree
    d = d_list(di);           % Get current degree
    b = b_func(d);            % Right endpoint for this degree
    x_eval = linspace(a, b, n_pts)';  % Query grid for this degree

    x_uni = mesh_uniform(a, b, d);  % Uniform mesh
    x_ch1 = mesh_cheby1(a, b, d);  % Cheby1 mesh
    x_ch2 = mesh_cheby2(a, b, d);  % Cheby2 mesh

    [ky_uni, k1_uni, Lam_uni, H_uni] = run_conditioning(x_uni, d, x_eval);  % Uniform conditioning
    [ky_ch1, k1_ch1, Lam_ch1, H_ch1] = run_conditioning(x_ch1, d, x_eval);  % Cheby1 conditioning
    [ky_ch2, k1_ch2, Lam_ch2, H_ch2] = run_conditioning(x_ch2, d, x_eval);  % Cheby2 conditioning

    Lam_uni_all(di) = Lam_uni;  % Store for trend plot
    Lam_ch1_all(di) = Lam_ch1;  % Store for trend plot

    fprintf('%-6d %-10s %-14.4e %-14.4e %-14.4e %-14.4e\n', d, 'Uniform', ...  % Uniform row
        Lam_uni, H_uni, max(ky_uni(isfinite(ky_uni))), max(k1_uni));
    fprintf('%-6d %-10s %-14.4e %-14.4e %-14.4e %-14.4e\n', d, 'Cheby1',  ...  % Cheby1 row
        Lam_ch1, H_ch1, max(ky_ch1(isfinite(ky_ch1))), max(k1_ch1));
    fprintf('%-6d %-10s %-14.4e %-14.4e %-14.4e %-14.4e\n', d, 'Cheby2',  ...  % Cheby2 row
        Lam_ch2, H_ch2, max(ky_ch2(isfinite(ky_ch2))), max(k1_ch2));
    fprintf('%s\n', repmat('-', 1, 76));  % Separator between degrees

    % ---- Subtask 3: accuracy for this degree ----
    f_exact = product_eval(double(x_eval), 1.0, (1:d)');  % Exact f2 values

    [p_b2_uni, p_inc_uni, p_dec_uni, p_lej_uni] = eval_all_forms(x_uni, d, x_eval, 'uniform');  % Uniform forms
    [p_b2_ch1, p_inc_ch1, p_dec_ch1, p_lej_ch1] = eval_all_forms(x_ch1, d, x_eval, 'cheby1');  % Cheby1 forms

    err_b2_uni_all(di)  = max(abs(double(p_b2_uni)  - f_exact));  % Bary2 uniform error
    err_b2_ch1_all(di)  = max(abs(double(p_b2_ch1)  - f_exact));  % Bary2 Cheby1 error

    % Higham-style plots only for largest degree (d=25, ~30 points)
    if d == d_list(end)

        % Stability bound for Bary Form 2 (Set 8 Section 5.4 / Higham 2004):
        %   (3n+4)*u*kappa_y + (3n+2)*u*kappa_1, multiplied by |f_exact|
        u = double(u_sp);
        stab_uni = ((3*d+4)*u*ky_uni + (3*d+2)*u*k1_uni) .* abs(f_exact);  % Bary2 bound uniform
        stab_ch1 = ((3*d+4)*u*ky_ch1 + (3*d+2)*u*k1_ch1) .* abs(f_exact);  % Bary2 bound Cheby1

        [p_b2_ch2, p_inc_ch2, p_dec_ch2, p_lej_ch2] = eval_all_forms(x_ch2, d, x_eval, 'cheby2');  % Cheby2 forms

        figure(3); clf; hold on; grid on;  % Figure 3: uniform Higham-style
        semilogy(x_eval, max(abs(double(p_b2_uni)  - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Bary2');       % Bary2
        semilogy(x_eval, max(abs(double(p_inc_uni) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton inc');  % Newton inc
        semilogy(x_eval, max(abs(double(p_dec_uni) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton dec');  % Newton dec
        semilogy(x_eval, max(abs(double(p_lej_uni) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton Leja'); % Newton Leja
        semilogy(x_eval, max(stab_uni, 1e-20), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Stability bound');              % Stability bound
        xlabel('x'); ylabel('|p(x) - f2(x)|');                                                  % Label axes
        title(sprintf('Figure 3: Pointwise error, f2, uniform mesh, d=%d', d));                 % Set title
        legend('Location','best'); hold off;                                                     % Add legend

        figure(4); clf; hold on; grid on;  % Figure 4: Cheby1 Higham-style
        semilogy(x_eval, max(abs(double(p_b2_ch1)  - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Bary2');       % Bary2
        semilogy(x_eval, max(abs(double(p_inc_ch1) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton inc');  % Newton inc
        semilogy(x_eval, max(abs(double(p_dec_ch1) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton dec');  % Newton dec
        semilogy(x_eval, max(abs(double(p_lej_ch1) - f_exact), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton Leja'); % Newton Leja
        semilogy(x_eval, max(stab_ch1, 1e-20), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Stability bound');              % Stability bound
        xlabel('x'); ylabel('|p(x) - f2(x)|');                                                  % Label axes
        title(sprintf('Figure 4: Pointwise error, f2, Cheby1 mesh, d=%d', d));                  % Set title
        legend('Location','best'); hold off;                                                     % Add legend

        figure(1); clf; hold on; grid on;  % Figure 1: Lebesgue function largest degree
        plot(x_eval, k1_uni, 'LineWidth', 1.5, 'DisplayName', sprintf('Uniform (\\Lambda_n=%.2e)', Lam_uni));  % Uniform
        plot(x_eval, k1_ch1, 'LineWidth', 1.5, 'DisplayName', sprintf('Cheby1  (\\Lambda_n=%.2e)', Lam_ch1));  % Cheby1
        plot(x_eval, k1_ch2, 'LineWidth', 1.5, 'DisplayName', sprintf('Cheby2  (\\Lambda_n=%.2e)', Lam_ch2));  % Cheby2
        xlabel('x'); ylabel('\kappa(x,n,1)');                                                    % Label axes
        title(sprintf('Figure 1: Lebesgue function, f2, d=%d', d));                             % Set title
        legend('Location','best'); hold off;                                                     % Add legend

        figure(2); clf; hold on; grid on;  % Figure 2: kappa(x,n,y) largest degree
        semilogy(x_eval, max(ky_uni, 1e-20), 'LineWidth', 1.5, 'DisplayName', sprintf('Uniform (H_n=%.2e)', H_uni));  % Uniform
        semilogy(x_eval, max(ky_ch1, 1e-20), 'LineWidth', 1.5, 'DisplayName', sprintf('Cheby1  (H_n=%.2e)', H_ch1));  % Cheby1
        semilogy(x_eval, max(ky_ch2, 1e-20), 'LineWidth', 1.5, 'DisplayName', sprintf('Cheby2  (H_n=%.2e)', H_ch2));  % Cheby2
        xlabel('x'); ylabel('\kappa(x,n,y)');                                                    % Label axes
        title(sprintf('Figure 2: kappa(x,n,y), f2, d=%d', d));                                  % Set title
        legend('Location','best'); hold off;                                                     % Add legend
    end
end
fprintf('\n');  % Blank line after table

% ============================================================
% SUBTASK 3: Summary error table across all degrees
% ============================================================

fprintf('------------------------------------------------------------\n');   % Print separator
fprintf('SUBTASK 3: Accuracy summary across all degrees\n');                  % Section header
fprintf('------------------------------------------------------------\n\n'); % Print separator

fprintf('Table 2: max |p(x) - f2(x)| for Bary2 (summary across degrees)\n');  % Table header
fprintf('%-6s %-18s %-18s\n', 'd', 'Bary2 Uniform', 'Bary2 Cheby1');           % Column headers
fprintf('%s\n', repmat('-', 1, 44));                                             % Separator
for di = 1 : length(d_list)  % Loop over each degree
    fprintf('%-6d %-18.4e %-18.4e\n', d_list(di), err_b2_uni_all(di), err_b2_ch1_all(di));  % Print row
end
fprintf('\n');  % Blank line

% Figure 5: Lambda_n vs degree
figure(5); clf; hold on; grid on;                                                                      % Create figure
plot(d_list, Lam_uni_all, '-o', 'LineWidth', 1.5, 'DisplayName', 'Uniform');                          % Uniform Lambda_n
plot(d_list, Lam_ch1_all, '-s', 'LineWidth', 1.5, 'DisplayName', 'Cheby1');                           % Cheby1 Lambda_n
xticks(d_list);                                                                                        % Set x ticks
xlabel('degree d'); ylabel('\Lambda_n');                                                               % Label axes
title('Figure 5: Lebesgue constant \Lambda_n vs degree, f2');                                         % Set title
legend('Location','best'); hold off;                                                                   % Add legend

% Figure 6: Max error vs degree for Bary2
figure(6); clf; hold on; grid on;                                                                      % Create figure
semilogy(d_list, err_b2_uni_all, '-o', 'LineWidth', 1.5, 'DisplayName', 'Bary2 Uniform');             % Uniform errors
semilogy(d_list, err_b2_ch1_all, '-s', 'LineWidth', 1.5, 'DisplayName', 'Bary2 Cheby1');             % Cheby1 errors
xticks(d_list);                                                                                        % Set x ticks
xlabel('degree d'); ylabel('max |p(x) - f2(x)|');                                                     % Label axes
title('Figure 6: Max error vs degree for Bary2, f2');                                                 % Set title
legend('Location','best'); hold off;                                                                   % Add legend

% ============================================================
% OVERDETERMINED TESTS
% ============================================================

fprintf('------------------------------------------------------------\n');   % Print separator
fprintf('OVERDETERMINED TESTS: more points than degree+1 for f2\n');         % Section header
fprintf('------------------------------------------------------------\n\n'); % Print separator

fprintf('Table 3: Overdetermined interpolation error (f2, true degree=%d)\n', d_over_deg);  % Table header
fprintf('%-8s %-10s %-20s %-20s\n', 'n_pts', 'Mesh', 'max |p-f2| Bary2', 'max |p-f2| Newton inc');  % Column headers
fprintf('%s\n', repmat('-', 1, 62));  % Separator

b_over  = b_func(d_over_deg);                       % Right endpoint for overdetermined tests
x_eval_over = linspace(a, b_over, n_pts)';          % Query grid for overdetermined tests
f_exact_over = product_eval(double(x_eval_over), 1.0, (1:d_over_deg)');  % Exact f2 values

for oi = 1 : length(d_over_list)  % Loop over each overdetermined size
    n_over = d_over_list(oi);      % Number of interpolation points minus 1

    x_uni_o = mesh_uniform(a, b_over, n_over);  % Uniform overdetermined mesh
    x_ch1_o = mesh_cheby1(a, b_over, n_over);   % Cheby1 overdetermined mesh

    for mi = 1 : 2  % Loop over mesh types
        if mi == 1
            x_m = x_uni_o;  mesh_name = 'Uniform';  mesh_type = 'uniform';  % Select uniform mesh
        else
            x_m = x_ch1_o;  mesh_name = 'Cheby1';   mesh_type = 'cheby1';   % Select Cheby1 mesh
        end

        x_s   = single(x_m);                                                 % Cast nodes to single
        fv    = f2(x_s, d_over_deg);                                         % f2 values at overdetermined nodes
        beta  = bary2_coeffs(double(x_m), n_over, mesh_type);                % Bary2 weights O(n)

        p_b2  = bary2_eval(x_s, fv, single(beta), single(x_eval_over));      % Bary2 evaluation
        c_inc = newton_coeffs(x_s, fv);                                      % Newton coeffs increasing
        p_inc = newton_eval(x_s, c_inc, single(x_eval_over));                % Newton increasing eval

        err_b2  = max(abs(double(p_b2)  - f_exact_over));  % Bary2 error
        err_inc = max(abs(double(p_inc) - f_exact_over));  % Newton inc error

        fprintf('%-8d %-10s %-20.4e %-20.4e\n', n_over+1, mesh_name, err_b2, err_inc);  % Print row
    end
end
fprintf('\n');  % Blank line

fprintf('DONE. Produced Tables 1-3 and Figures 1-6.\n');  % Print completion message
