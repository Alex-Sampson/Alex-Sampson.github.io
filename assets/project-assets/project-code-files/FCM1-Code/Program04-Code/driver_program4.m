function driver_program4()
% DRIVER_PROGRAM4  Main driver for Program 4 (Householder-QR least squares).
%
% Sections:
%   1) Structured least-squares tests via gen_structured_lls
%      - compare lsq_householder vs lsq_householder_incremental
%   2) Growth matrix: behavior of entries in the R factor vs n
%   3) Arithmetic mean as a 1D least-squares problem
%   4) Polynomial regression:
%        - monomial vs Chebyshev bases
%        - multiple point sets
%        - Gram matrices B'*B
%
% This driver assumes the following routines are in the same file:
%   lsq_householder.m
%   lsq_householder_incremental.m
%   my_qr_householder.m
%   apply_Qt_to_vec.m
%   hh_vec.m
%   gen_structured_lls.m
%   growth_matrix.m
%   chebyshev_nodes.m
%   build_poly_matrix.m
%   gram_matrix.m

    rng(1);  % use the same random seed to reproduce randomness

    fprintf('\n==============================================\n');
    fprintf(' Program 4 Driver: Householder QR Least Squares\n');
    fprintf('==============================================\n\n');

    % Toggle plotting on/off
    make_plots = true;

    % 1) Structured least-squares tests
    run_structured_lls_tests(make_plots);

    % 2) Growth matrix and R-factor behavior
    run_growth_matrix_tests(make_plots);

    % 3) Arithmetic mean as least squares
    run_arithmetic_mean_tests(make_plots);

    % 4) Polynomial regression experiments
    run_poly_regression_tests(make_plots);

    fprintf('\nDriver finished.\n');
end


% ======================================================================
% 1) Structured LS tests (cleaner plots, split by kind)
% ======================================================================

function run_structured_lls_tests(make_plots)

    fprintf('----------------------------------------------\n');
    fprintf('1) Structured least-squares tests\n');
    fprintf('   (gen_structured_lls + lsq_householder vs incremental)\n');
    fprintf('----------------------------------------------\n');

    % kind = 1 : square, nonsingular (n = k)
    % kind = 2 : tall, consistent
    % kind = 3 : tall, inconsistent
    kinds   = [1 2 3];

    % Problem size range
    k_list  = 10:10:200;

    dscale  = 1e-1;   % moderate inconsistency
    nref    = 6;      % number of random reflectors in gen_structured_lls
    tol     = 1e-12;

    fprintf('  %-6s %-4s %-4s %-12s %-12s %-12s %-12s\n', ...
        'kind', 'n', 'k', 'rel||x1-x2||', 'rel||c1-c2||', ...
        'rel||d1-d2||', 'rel resid (x1)');

    % collect data per kind for plots
    n_list_kind   = cell(3,1);
    relx_kind     = cell(3,1);
    relres_kind   = cell(3,1);

    for k = k_list
        for kind = kinds
            if kind == 1
                n = k;
            else
                n = 2*k;
            end

            [A, b, xtrue, Rtrue, ctrue, dtrue] = gen_structured_lls(n, k, kind, dscale, nref);

            [x1, R1, c1, d1] = lsq_householder(A, b, tol); 
            [x2, R2, c2, d2] = lsq_householder_incremental(A, b, tol); 

            rel_x = rel_diff(x1, x2);
            rel_c = rel_diff(c1, c2);
            rel_d = rel_diff(d1, d2);

            rel_resid_x1 = norm(b - A*x1, 2) / max(1, norm(b, 2));

            fprintf('  %-6d %-4d %-4d %-12.3e %-12.3e %-12.3e %-12.3e\n', ...
                kind, n, k, rel_x, rel_c, rel_d, rel_resid_x1);

            % store for this kind
            idx = kind;
            n_list_kind{idx}(end+1,1) = n;
            relx_kind{idx}(end+1,1)   = rel_x;
            relres_kind{idx}(end+1,1) = rel_resid_x1;
        end
    end

    fprintf('\n(Expect all lsq_householder vs incremental differences to be at or near roundoff.)\n\n');

    if make_plots
        markers = {'o-','s-','^-'};  % one style per kind

%------------- Plot rel||x1-x2|| vs n, one curve per kind --------------------

        figure; hold on;
        for kind = kinds
            n_vals  = n_list_kind{kind};
            y_vals  = relx_kind{kind};
            y_vals(y_vals < eps) = eps;   % clip insanely small values so they don't underflow on log scale
            loglog(n_vals, y_vals, markers{kind}, 'LineWidth', 1.0);
        end

        grid on;
        xlabel('n');
        ylabel('rel ||x_{full} - x_{inc}||_2');
        title('Program 4: Difference between full and incremental LS solutions');
        legend('kind = 1 (square)', 'kind = 2 (tall, cons.)', ...
               'kind = 3 (tall, incons.)', 'Location', 'southwest');
        hold off;

%-------------- Plot relative residual vs n, one curve per kind--------------

        figure; hold on;
        for kind = kinds
            n_vals = n_list_kind{kind};
            y_vals = relres_kind{kind};
            loglog(n_vals, y_vals, markers{kind}, 'LineWidth', 1.0);
        end
        grid on;
        xlabel('n');
        ylabel('rel residual for x_{full}');
        title('Program 4: Relative residual vs n for structured LS problems');
        legend('kind = 1 (square)', 'kind = 2 (tall, cons.)', ...
               'kind = 3 (tall, incons.)', 'Location', 'northwest');
        hold off;
    end
end


% ======================================================================
% 2) Growth matrix and R-factor behavior
% ======================================================================

function run_growth_matrix_tests(make_plots)
    fprintf('----------------------------------------------\n');
    fprintf('2) Growth matrix: R entries from Householder QR\n');
    fprintf('----------------------------------------------\n');

    % Growth matrix from Program 3:
    %   A(i,i) = 1
    %   A(i,j) = -1 for j < i
    %   A(i,n) = 1
    %
    % Here we factor it with Householder QR and look at the size of entries in R.

    % n from 10 up to 500 in steps of 10

    n_list = 10:10:500;

    maxR   = zeros(length(n_list), 1);
    normR  = zeros(length(n_list), 1);

    fprintf('  %-6s %-12s %-12s\n', 'n', 'max|R_ij|', '||R||_2 (cond. approx)');

    for idx = 1:length(n_list)
        n = n_list(idx);
        A = growth_matrix(n);

        [W, tau] = my_qr_householder(A); 
        R = triu(W(1:n, 1:n));

        maxR(idx)  = max(abs(R(:)));
        normR(idx) = norm(R, 2);

        fprintf('  %-6d %-12.3e %-12.3e\n', n, maxR(idx), normR(idx));
    end

    if make_plots
        figure;
        semilogy(n_list, maxR, 'o-','LineWidth',1.0);
        grid on;
        xlabel('n');
        ylabel('max |R_{ij}| (growth matrix)');
        title('Program 4: Magnitude of entries in R for the growth matrix');
    end

    fprintf('\n(Use these values to compare with the U factors from LU in Program 3.)\n\n');
end


% ======================================================================
% 3) Arithmetic mean as least squares (cleaner plot)
% ======================================================================

function run_arithmetic_mean_tests(make_plots)
    fprintf('----------------------------------------------\n');
    fprintf('3) Arithmetic mean as a least-squares problem\n');
    fprintf('----------------------------------------------\n');

    % min_beta ||eta - 1*beta||_2  with A = ones(n,1)

    % sample sizes
    n_list = [10 20 50 100 200 500 1000 2000 5000];

    tol    = 1e-12;

    fprintf('  %-6s %-12s %-12s %-12s\n', 'n', 'beta_true', 'beta_LS', '|beta_LS-beta_true|');

    err_all = zeros(length(n_list),1);

    for idx = 1:length(n_list)
        n   = n_list(idx);
        eta = randn(n, 1);
        A   = ones(n, 1);
        b   = eta;

        beta_true = mean(eta);

        [x_hat, ~, c, d] = lsq_householder(A, b, tol); 

        beta_ls   = x_hat;
        err_beta  = abs(beta_ls - beta_true);
        err_all(idx) = err_beta;

        fprintf('  %-6d %-12.5f %-12.5f %-12.3e\n', n, beta_true, beta_ls, err_beta);
    end

    fprintf('\n(Arithmetic mean is recovered as the 1D LS solution up to rounding.)\n\n');

    if make_plots
        figure;
        semilogx(n_list, err_all, 'o-','LineWidth',1.0);
        grid on;
        xlabel('n');
        ylabel('|beta_{LS} - beta_{true}|');
        title('Program 4: Error in LS-based arithmetic mean vs n');
    end
end


% ======================================================================
% 4) Polynomial regression: monomial vs Chebyshev
% ======================================================================

function run_poly_regression_tests(make_plots)
    fprintf('----------------------------------------------\n');
    fprintf('4) Polynomial regression (monomial vs Chebyshev bases)\n');
    fprintf('----------------------------------------------\n');

    degrees       = [1 2 3];
    tol           = 1e-12;

    % Larger uniform meshes
    n_uniform1    = 50;
    n_uniform2    = 100;

    % Fine grid for comparing polynomials
    x_fine = linspace(-1, 1, 200)';

    fprintf('  %-4s %-10s %-8s %-14s %-14s %-14s\n', ...
        'd', 'mesh', 'n', 'cond(G_mono)', 'cond(G_cheb)', 'max|p_m-p_c|');

    for d = degrees
        a_true_full = [1; -0.5; 0.25; 0.1];
        a_true      = a_true_full(1:(d+1));

        % --- Mesh 1: uniform points, n_uniform1 ---
        x1 = linspace(-1, 1, n_uniform1)';
        run_poly_case(d, x1, a_true, x_fine, tol, 'uniform-50');

        % --- Mesh 2: uniform points, n_uniform2 ---
        x2 = linspace(-1, 1, n_uniform2)';
        run_poly_case(d, x2, a_true, x_fine, tol, 'uniform-100');

        % --- Mesh 3: Chebyshev nodes (d+1 points) ---
        x3 = chebyshev_nodes(d);  % d+1 Chebyshev nodes in [-1,1]
        run_poly_case(d, x3, a_true, x_fine, tol, 'Cheb-nodes');
    end

    if make_plots
        % Example plot: for one specific case, visualize the fits.

        d = 3;
        x_plot = linspace(-1, 1, 100)';

        a_true_full = [1; -0.5; 0.25; 0.1];
        a_true      = a_true_full(1:(d+1));

        x_case = linspace(-1, 1, 100)';   % match the larger uniform mesh
        Bm     = build_poly_matrix(x_case, d, 0);
        y_case = Bm * a_true;

        a_mono = lsq_householder(Bm, y_case, tol);
        Bc     = build_poly_matrix(x_case, d, 1);
        a_cheb = lsq_householder(Bc, y_case, tol);

        Bm_plot = build_poly_matrix(x_plot, d, 0);
        Bc_plot = build_poly_matrix(x_plot, d, 1);
        p_true  = Bm_plot * a_true;
        p_mono  = Bm_plot * a_mono;
        p_cheb  = Bc_plot * a_cheb;

        figure;
        plot(x_plot, p_true, 'k-', 'LineWidth',1.5); hold on;
        plot(x_plot, p_mono, 'b--', 'LineWidth',1.2);
        plot(x_plot, p_cheb, 'r-.', 'LineWidth',1.2);
        grid on;
        xlabel('x'); ylabel('p(x)');
        title(sprintf('Degree d=%d regression: monomial vs Chebyshev (n=100 uniform)', d));
        legend('true', 'monomial LS', 'Chebyshev LS', 'Location','best');
        hold off;
    end

    fprintf('\n(Polynomials from monomial vs Chebyshev bases should agree up to rounding.\n');
    fprintf('  Compare cond(G) across meshes/bases in your report.)\n\n');
end


% ======================================================================
% one polynomial regression case
% ======================================================================

function run_poly_case(d, x, a_true, x_fine, tol, mesh_label)
    n = length(x);

    B_mono = build_poly_matrix(x, d, 0);
    y      = B_mono * a_true;

    B_cheb = build_poly_matrix(x, d, 1);

    a_mono = lsq_householder(B_mono, y, tol);
    a_cheb = lsq_householder(B_cheb, y, tol);

    G_mono = gram_matrix(B_mono);
    G_cheb = gram_matrix(B_cheb);

    cond_Gm = cond(G_mono);
    cond_Gc = cond(G_cheb);

    Bm_fine = build_poly_matrix(x_fine, d, 0);
    Bc_fine = build_poly_matrix(x_fine, d, 1);
    p_mono  = Bm_fine * a_mono;
    p_cheb  = Bc_fine * a_cheb;

    max_diff = max(abs(p_mono - p_cheb));

    fprintf('  %-4d %-10s %-8d %-14.3e %-14.3e %-14.3e\n', ...
        d, mesh_label, n, cond_Gm, cond_Gc, max_diff);
end


% ======================================================================
% Small utility: relative difference between two vectors
% ======================================================================

function r = rel_diff(u, v)
    if isempty(u) && isempty(v)
        r = 0;
        return;
    end
    if isempty(u) || isempty(v)
        r = Inf;
        return;
    end
    u = u(:); v = v(:);
    denom = max(1, norm(v, 2));
    r = norm(u - v, 2) / denom;
end
