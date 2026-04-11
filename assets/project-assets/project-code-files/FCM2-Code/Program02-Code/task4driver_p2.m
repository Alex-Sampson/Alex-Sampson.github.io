% task4driver_p2.m
%
% Task 4: Interpolation of f3(x) = l_n(x), the n-th Lagrange basis function.
%
%   f3(x_i) = 0,  0 <= i <= n-1
%   f3(x_n) = 1
%
% Key distinction from Tasks 2 & 3: no product form is available for f3
% since f3 depends on the mesh itself.  The "exact" denominator for
% condition numbers is therefore taken from bary1_eval_dp (double precision).
%
% Subroutines used:
%   mesh_uniform, mesh_cheby1, mesh_cheby2, leja_order
%   bary1_coeffs, bary1_eval_dp
%   bary2_eval, newton_coeffs, newton_eval
%   cond_kappa, lebesgue_lambda, hilbert_H
%   f3
%
% Subtask 1: Interpolation on uniform, Cheby1, Cheby2 meshes for n=21,25,29
%            including overdetermined tests
% Subtask 2: Conditioning - kappa(x,n,y), kappa(x,n,1), Lambda_n, H_n
% Subtask 3: Accuracy and stability of single precision forms vs exact
%            Higham-style plots for n=29 on uniform and Cheby1 meshes
%
% Output:
%   Table 1: Lambda_n and H_n for all degrees and mesh types
%   Table 2: Max error for all forms, all degrees, all mesh types
%   Table 3: Overdetermined test errors
%   Figure 1: Lebesgue function, uniform vs Cheby1, n=29
%   Figure 2: kappa(x,n,y), uniform vs Cheby1, n=29
%   Figure 3: Pointwise error all forms, uniform mesh, n=29 (Higham-style)
%   Figure 4: Pointwise error all forms, Cheby1 mesh, n=29 (Higham-style)
%   Figure 5: Lambda_n vs n for uniform and Cheby1
%   Figure 6: Max error vs n for Bary2 on uniform and Cheby1

clear; clc; close all;  % Clear workspace, command window, and close all figures

% ----------------  settings  ----------------
n_list      = [21, 25, 29];   % Degrees to test (>=2 above 20, n=29 required for Higham)
n_over_list = [35, 45];       % Overdetermined: more points than degree+1
n_over_deg  = 21;             % True degree of f3 for overdetermined tests
a           = -1.0;           % Left endpoint (standard interval for Lagrange basis)
b           =  1.0;           % Right endpoint
n_pts       = 1000;           % Number of query points
u_sp        = single(eps('single'));  % Unit roundoff for single precision

x_eval = linspace(a, b, n_pts)';  % Dense query grid over [a, b]

fprintf('\n============================================================\n');  % Print separator
fprintf('TASK 4: f3(x) = l_n(x), n-th Lagrange basis function\n');           % Print task header
fprintf('Interval [%.4g, %.4g],  degrees: ', a, b);                          % Print interval
fprintf('%d ', n_list);  fprintf('\n');                                        % Print degree list
fprintf('Note: no product form available; bary1_eval_dp used as reference\n');% Print key note
fprintf('============================================================\n\n');  % Print separator

% ============================================================
% LOCAL HELPERS
% ============================================================

function [p_b2, p_inc, p_dec, p_lej] = eval_all_forms(x_nodes, x_eval, mesh_type)
    % Evaluate all four single-precision forms for f3 on a given mesh.
    n      = length(x_nodes) - 1;                            % Degree
    x_s    = single(x_nodes);                                % Cast nodes to single
    fv     = f3(x_s, x_s);                                   % f3 values at nodes (single)
    beta   = bary2_coeffs(double(x_nodes), n, mesh_type);    % Bary2 weights O(n)

    p_b2   = bary2_eval(x_s, fv, single(beta), single(x_eval));  % Barycentric Form 2

    c_inc  = newton_coeffs(x_s, fv);                               % Newton coeffs increasing
    p_inc  = newton_eval(x_s, c_inc, single(x_eval));              % Newton increasing

    x_dec  = flipud(x_s);  fv_dec = flipud(fv);                    % Reverse order
    c_dec  = newton_coeffs(x_dec, fv_dec);                         % Newton coeffs decreasing
    p_dec  = newton_eval(x_dec, c_dec, single(x_eval));            % Newton decreasing

    x_lej  = single(leja_order(x_nodes));                          % Leja ordering
    fv_lej = f3(x_lej, x_s);                                       % f3 values at Leja nodes
    c_lej  = newton_coeffs(x_lej, fv_lej);                         % Newton coeffs Leja
    p_lej  = newton_eval(x_lej, c_lej, single(x_eval));            % Newton Leja
end

function [kappa_y, kappa_1, Lambda_n, H_n, f_exact] = run_conditioning(x_nodes, x_eval)
    % Compute condition numbers and exact reference for f3 on a given mesh.
    % No product form: bary1_eval_dp provides both p_exact and the denominator.
    n      = length(x_nodes) - 1;                                   % Degree
    fv_dp  = double(f3(double(x_nodes), double(x_nodes)));          % f3 values in double
    gamma  = bary1_coeffs(double(x_nodes), n, ...
                          @(t) double(f3(t, double(x_nodes))));      % Bary1 weights (double)

    [f_exact, ~, denom] = bary1_eval_dp(double(x_nodes), fv_dp, ...
                                         gamma, double(x_eval));     % Exact values and denominator

    [kappa_y, kappa_1]  = cond_kappa(double(x_nodes), fv_dp, gamma, ...
                                      double(x_eval), abs(denom));   % Condition numbers

    Lambda_n = lebesgue_lambda(kappa_1);  % Lebesgue constant
    H_n      = hilbert_H(kappa_y);        % H_n norm
end

% ============================================================
% SUBTASKS 1 & 2: Multi-degree conditioning sweep
% ============================================================

fprintf('------------------------------------------------------------\n');   % Print separator
fprintf('SUBTASKS 1 & 2: Conditioning sweep over degrees\n');                % Section header
fprintf('------------------------------------------------------------\n\n'); % Print separator

fprintf('Table 1: Conditioning summary for f3 (all degrees and meshes)\n');             % Table header
fprintf('%-6s %-10s %-14s %-14s %-14s %-14s\n', ...                                    % Column headers
    'n', 'Mesh', 'Lambda_n', 'H_n', 'max kappa_y', 'max kappa_1');
fprintf('%s\n', repmat('-', 1, 76));  % Separator

Lam_uni_all     = zeros(length(n_list), 1);  % Lambda_n uniform storage
Lam_ch1_all     = zeros(length(n_list), 1);  % Lambda_n Cheby1 storage
err_b2_uni_all  = zeros(length(n_list), 1);  % Bary2 uniform error storage
err_b2_ch1_all  = zeros(length(n_list), 1);  % Bary2 Cheby1 error storage

for ni = 1 : length(n_list)  % Loop over each degree
    n = n_list(ni);           % Get current degree

    x_uni = mesh_uniform(a, b, n);  % Uniform mesh
    x_ch1 = mesh_cheby1(a, b, n);  % Cheby1 mesh
    x_ch2 = mesh_cheby2(a, b, n);  % Cheby2 mesh

    [ky_uni, k1_uni, Lam_uni, H_uni, fex_uni] = run_conditioning(x_uni, x_eval);  % Uniform
    [ky_ch1, k1_ch1, Lam_ch1, H_ch1, fex_ch1] = run_conditioning(x_ch1, x_eval);  % Cheby1
    [ky_ch2, k1_ch2, Lam_ch2, H_ch2, ~      ] = run_conditioning(x_ch2, x_eval);  % Cheby2

    Lam_uni_all(ni) = Lam_uni;  % Store for trend plot
    Lam_ch1_all(ni) = Lam_ch1;  % Store for trend plot

    fprintf('%-6d %-10s %-14.4e %-14.4e %-14.4e %-14.4e\n', n, 'Uniform', ...  % Uniform row
        Lam_uni, H_uni, max(ky_uni(isfinite(ky_uni))), max(k1_uni));
    fprintf('%-6d %-10s %-14.4e %-14.4e %-14.4e %-14.4e\n', n, 'Cheby1',  ...  % Cheby1 row
        Lam_ch1, H_ch1, max(ky_ch1(isfinite(ky_ch1))), max(k1_ch1));
    fprintf('%-6d %-10s %-14.4e %-14.4e %-14.4e %-14.4e\n', n, 'Cheby2',  ...  % Cheby2 row
        Lam_ch2, H_ch2, max(ky_ch2(isfinite(ky_ch2))), max(k1_ch2));
    fprintf('%s\n', repmat('-', 1, 76));  % Separator between degrees

    % ---- Subtask 3: accuracy for this degree ----
    [p_b2_uni, p_inc_uni, p_dec_uni, p_lej_uni] = eval_all_forms(x_uni, x_eval, 'uniform');  % Uniform forms
    [p_b2_ch1, p_inc_ch1, p_dec_ch1, p_lej_ch1] = eval_all_forms(x_ch1, x_eval, 'cheby1');  % Cheby1 forms

    err_b2_uni_all(ni) = max(abs(double(p_b2_uni) - fex_uni));  % Bary2 uniform error
    err_b2_ch1_all(ni) = max(abs(double(p_b2_ch1) - fex_ch1));  % Bary2 Cheby1 error

    % Higham-style plots only for n=29 (30 points, matches Higham paper)
    if n == 29

        % Stability bound for Bary Form 2 (Set 8 Section 5.4 / Higham 2004):
        %   (3n+4)*u*kappa_y + (3n+2)*u*kappa_1, multiplied by |f_exact|
        u = double(u_sp);
        stab_uni = ((3*n+4)*u*ky_uni + (3*n+2)*u*k1_uni) .* abs(fex_uni);  % Bary2 bound uniform
        stab_ch1 = ((3*n+4)*u*ky_ch1 + (3*n+2)*u*k1_ch1) .* abs(fex_ch1);  % Bary2 bound Cheby1

        [p_b2_ch2, p_inc_ch2, p_dec_ch2, p_lej_ch2] = eval_all_forms(x_ch2, x_eval, 'cheby2');  % Cheby2 forms
        [~, ~, ~, ~, fex_ch2] = run_conditioning(x_ch2, x_eval);  % Cheby2 exact reference

        figure(1); clf; hold on; grid on;  % Figure 1: Lebesgue function n=29
        plot(x_eval, k1_uni, 'LineWidth', 1.5, 'DisplayName', sprintf('Uniform (\\Lambda_n=%.2e)', Lam_uni));  % Uniform
        plot(x_eval, k1_ch1, 'LineWidth', 1.5, 'DisplayName', sprintf('Cheby1  (\\Lambda_n=%.2e)', Lam_ch1));  % Cheby1
        plot(x_eval, k1_ch2, 'LineWidth', 1.5, 'DisplayName', sprintf('Cheby2  (\\Lambda_n=%.2e)', Lam_ch2));  % Cheby2
        xlabel('x'); ylabel('\kappa(x,n,1)');                                                  % Label axes
        title(sprintf('Figure 1: Lebesgue function, f3, n=%d', n));                           % Set title
        legend('Location','best'); hold off;                                                   % Add legend

        figure(2); clf; hold on; grid on;  % Figure 2: kappa(x,n,y) n=29
        semilogy(x_eval, max(ky_uni, 1e-20), 'LineWidth', 1.5, 'DisplayName', sprintf('Uniform (H_n=%.2e)', H_uni));  % Uniform
        semilogy(x_eval, max(ky_ch1, 1e-20), 'LineWidth', 1.5, 'DisplayName', sprintf('Cheby1  (H_n=%.2e)', H_ch1));  % Cheby1
        semilogy(x_eval, max(ky_ch2, 1e-20), 'LineWidth', 1.5, 'DisplayName', sprintf('Cheby2  (H_n=%.2e)', H_ch2));  % Cheby2
        xlabel('x'); ylabel('\kappa(x,n,y)');                                                  % Label axes
        title(sprintf('Figure 2: kappa(x,n,y), f3, n=%d', n));                                % Set title
        legend('Location','best'); hold off;                                                   % Add legend

        figure(3); clf; hold on; grid on;  % Figure 3: Higham-style uniform n=29
        semilogy(x_eval, max(abs(double(p_b2_uni)  - fex_uni), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Bary2');       % Bary2
        semilogy(x_eval, max(abs(double(p_inc_uni) - fex_uni), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton inc');  % Newton inc
        semilogy(x_eval, max(abs(double(p_dec_uni) - fex_uni), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton dec');  % Newton dec
        semilogy(x_eval, max(abs(double(p_lej_uni) - fex_uni), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton Leja'); % Newton Leja
        semilogy(x_eval, max(stab_uni, 1e-20), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Stability bound');              % Stability bound
        xlabel('x'); ylabel('|p(x) - f3(x)|');                                                % Label axes
        title(sprintf('Figure 3: Pointwise error, f3, uniform mesh, n=%d', n));               % Set title
        legend('Location','best'); hold off;                                                   % Add legend

        figure(4); clf; hold on; grid on;  % Figure 4: Higham-style Cheby1 n=29
        semilogy(x_eval, max(abs(double(p_b2_ch1)  - fex_ch1), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Bary2');       % Bary2
        semilogy(x_eval, max(abs(double(p_inc_ch1) - fex_ch1), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton inc');  % Newton inc
        semilogy(x_eval, max(abs(double(p_dec_ch1) - fex_ch1), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton dec');  % Newton dec
        semilogy(x_eval, max(abs(double(p_lej_ch1) - fex_ch1), 1e-20), 'LineWidth', 1.5, 'DisplayName', 'Newton Leja'); % Newton Leja
        semilogy(x_eval, max(stab_ch1, 1e-20), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Stability bound');              % Stability bound
        xlabel('x'); ylabel('|p(x) - f3(x)|');                                                % Label axes
        title(sprintf('Figure 4: Pointwise error, f3, Cheby1 mesh, n=%d', n));                % Set title
        legend('Location','best'); hold off;                                                   % Add legend
    end
end
fprintf('\n');  % Blank line after table

% ============================================================
% SUBTASK 3: Summary error table
% ============================================================

fprintf('------------------------------------------------------------\n');   % Print separator
fprintf('SUBTASK 3: Accuracy summary across all degrees\n');                  % Section header
fprintf('------------------------------------------------------------\n\n'); % Print separator

fprintf('Table 2: max |p(x) - f3(x)| for Bary2 (summary across degrees)\n');  % Table header
fprintf('%-6s %-18s %-18s\n', 'n', 'Bary2 Uniform', 'Bary2 Cheby1');           % Column headers
fprintf('%s\n', repmat('-', 1, 44));                                             % Separator
for ni = 1 : length(n_list)  % Loop over each degree
    fprintf('%-6d %-18.4e %-18.4e\n', n_list(ni), err_b2_uni_all(ni), err_b2_ch1_all(ni));  % Print row
end
fprintf('\n');  % Blank line

% Figure 5: Lambda_n vs degree
figure(5); clf; hold on; grid on;                                                                   % Create figure
plot(n_list, Lam_uni_all, '-o', 'LineWidth', 1.5, 'DisplayName', 'Uniform');                       % Uniform Lambda_n
plot(n_list, Lam_ch1_all, '-s', 'LineWidth', 1.5, 'DisplayName', 'Cheby1');                        % Cheby1 Lambda_n
xticks(n_list);                                                                                     % Set x ticks
xlabel('degree n'); ylabel('\Lambda_n');                                                            % Label axes
title('Figure 5: Lebesgue constant \Lambda_n vs degree, f3');                                      % Set title
legend('Location','best'); hold off;                                                                % Add legend

% Figure 6: Max error vs degree for Bary2
figure(6); clf; hold on; grid on;                                                                   % Create figure
semilogy(n_list, err_b2_uni_all, '-o', 'LineWidth', 1.5, 'DisplayName', 'Bary2 Uniform');          % Uniform errors
semilogy(n_list, err_b2_ch1_all, '-s', 'LineWidth', 1.5, 'DisplayName', 'Bary2 Cheby1');           % Cheby1 errors
xticks(n_list);                                                                                     % Set x ticks
xlabel('degree n'); ylabel('max |p(x) - f3(x)|');                                                  % Label axes
title('Figure 6: Max error vs degree for Bary2, f3');                                              % Set title
legend('Location','best'); hold off;                                                                % Add legend

% ============================================================
% OVERDETERMINED TESTS
% ============================================================

fprintf('------------------------------------------------------------\n');   % Print separator
fprintf('OVERDETERMINED TESTS: more points than degree+1 for f3\n');         % Section header
fprintf('------------------------------------------------------------\n\n'); % Print separator

x_uni_base = mesh_uniform(a, b, n_over_deg);  % Base mesh defining f3

fprintf('Table 3: Overdetermined interpolation error (f3, true degree=%d)\n', n_over_deg);  % Table header
fprintf('%-8s %-10s %-20s %-20s\n', 'n_pts', 'Mesh', 'max |p-f3| Bary2', 'max |p-f3| Newton inc');  % Column headers
fprintf('%s\n', repmat('-', 1, 62));  % Separator

[~, ~, ~, ~, f_exact_over] = run_conditioning(x_uni_base, x_eval);  % Exact f3 reference

for oi = 1 : length(n_over_list)  % Loop over overdetermined sizes
    n_over = n_over_list(oi);      % Number of interpolation points minus 1

    x_uni_o = mesh_uniform(a, b, n_over);  % Uniform overdetermined mesh
    x_ch1_o = mesh_cheby1(a, b, n_over);   % Cheby1 overdetermined mesh

    for mi = 1 : 2  % Loop over mesh types
        if mi == 1
            x_m = x_uni_o;  mesh_name = 'Uniform';  mesh_type = 'uniform';  % Select uniform mesh
        else
            x_m = x_ch1_o;  mesh_name = 'Cheby1';   mesh_type = 'cheby1';   % Select Cheby1 mesh
        end

        x_s    = single(x_m);                                                 % Cast to single
        fv     = f3(x_s, single(x_uni_base));                                 % f3 values at overdetermined nodes
        beta   = bary2_coeffs(double(x_m), n_over, mesh_type);                % Bary2 weights O(n)

        p_b2   = bary2_eval(x_s, fv, single(beta), single(x_eval));           % Bary2 evaluation
        c_inc  = newton_coeffs(x_s, fv);                                       % Newton coefficients
        p_inc  = newton_eval(x_s, c_inc, single(x_eval));                     % Newton evaluation

        err_b2  = max(abs(double(p_b2)  - f_exact_over));  % Bary2 error
        err_inc = max(abs(double(p_inc) - f_exact_over));  % Newton inc error

        fprintf('%-8d %-10s %-20.4e %-20.4e\n', n_over+1, mesh_name, err_b2, err_inc);  % Print row
    end
end
fprintf('\n');  % Blank line

fprintf('DONE. Produced Tables 1-3 and Figures 1-6.\n');  % Print completion message
