% driver_p4.m
%
% Program 4 Driver -- Composite Numerical Quadrature
%
% Works through all four tasks on all five test integrals using:
%   comp_trap          -- Composite Trapezoidal Rule (a priori sizing)
%   comp_trap_adaptive -- Composite Trapezoidal Rule (adaptive step-halving)
%   comp_gauss2        -- Composite 2-point Gauss-Legendre Rule
%
% Task i:   A priori step size prediction for CTR across a range of tolerances
% Task ii:  CTR with a priori Hm -- compare predicted vs observed vs true error
% Task iii: Adaptive CTR -- same tolerances, compare accuracy per eval
% Task iv:  CGL cost comparison -- evals to match CTR accuracy at O(H^4)
%
% Outputs:
%   Table 1:  Predicted Hm, m, and error bound per tolerance (Task i)
%   Table 2:  True error vs predicted bound -- tightness ratio (Task ii)
%   Table 3:  Adaptive vs a priori eval counts and eff. ratios (Task iii)
%   Table 4:  Three-way eval comparison CTR apri / CTR adap / CGL (Task iv)
%   Figure 1: Predicted vs true error vs m -- tightness of the CTR bound
%   Figure 2: Adaptive eval savings ratio vs tolerance (when adaptive pays off)
%   Figure 3: Convergence rates -- O(H^2) CTR vs O(H^4) CGL with reference slopes
%   Figure 4: Error vs total function evaluations -- efficiency per integrand
%   Figure 5: CGL vs CTR eval ratio vs tolerance -- where does CGL dominate?

clear; clc; close all;

% ---- color palette (clr_ prefix avoids clashes with output variable names) ----
clr_apri   = [0.122, 0.471, 0.706];   % steel blue   -- CTR a priori
clr_adap   = [0.890, 0.102, 0.110];   % crimson      -- CTR adaptive
clr_cgl    = [0.200, 0.627, 0.173];   % forest green -- CGL
clr_bound  = [0.596, 0.306, 0.639];   % purple       -- predicted bound
clr_true   = [0.000, 0.000, 0.000];   % black        -- true error
clr_ref2   = [0.5 0.5 0.5];           % grey         -- O(H^2) reference
clr_ref4   = [0.3 0.3 0.3];           % dark grey    -- O(H^4) reference

% ============================================================
% SETUP: test integrals, exact values, second derivative bounds
% ============================================================

% -- the five integrands -- (define the functions to integrate)
f1 = @(x) exp(x);                                % f1(x) = e^x
f2 = @(x) exp(sin(2*x)) .* cos(2*x);            % f2(x) = e^{sin(2x)} * cos(2x)
f3 = @(x) tanh(x);                               % f3(x) = tanh(x)
f4 = @(x) x .* cos(2*pi*x);                      % f4(x) = x * cos(2πx)
f5 = @(x) x + 1./x;                              % f5(x) = x + 1/x

funcs  = {f1, f2, f3, f4, f5};                   % put all functions into a cell array
fnames = {'f_1(x) = e^x', 'f_2(x) = e^{sin(2x)}cos(2x)', ...
          'f_3(x) = tanh(x)', 'f_4(x) = x cos(2\pi x)', 'f_5(x) = x + 1/x'}; % readable names
fshort = {'f1', 'f2', 'f3', 'f4', 'f5'};         % short labels for plots or filenames

% -- integration limits -- (left and right endpoints for each integral)
limits = [ 0,    3;
           0,    pi/3;
          -2,    1;
           0,    3.5;
           0.1,  2.5 ];

% -- exact values from the problem statement -- (precomputed true integrals)
exact(1) = exp(3) - 1;                           % exact integral of e^x from 0 to 3
exact(2) = 0.5 * (-1 + exp(sqrt(3)/2));         % exact value given for f2
exact(3) = log(cosh(1) / cosh(2));               % exact value given for f3
exact(4) = -1 / (2 * pi^2);                      % exact value given for f4
exact(5) = (2.5^2 - 0.1^2) / 2 + log(2.5 / 0.1);  % exact value for f5 over [0.1,2.5]

% -- M2 = max|f''(x)| estimated numerically on a dense grid --
% used only for a priori sizing -- fine-grained centered differences are
% plenty accurate for a bound that will be rounded up anyway
n_dense = 10000;                                 % number of points to sample for derivative est.
M2      = zeros(5, 1);                           % preallocate storage for max second derivative
for j = 1 : 5
    a_j = limits(j, 1);                          % left endpoint for j-th integral
    b_j = limits(j, 2);                          % right endpoint for j-th integral
    x_d = linspace(a_j, b_j, n_dense);           % dense grid across the interval
    h_d = (b_j - a_j) / (n_dense - 1);           % spacing between dense grid points
    fv  = arrayfun(funcs{j}, x_d);               % function values at the dense grid
    f2v = zeros(1, n_dense);                     % space for second derivative estimates
    % forward/backward difference at the endpoints (to avoid indexing issues)
    f2v(1)       = (fv(3) - 2*fv(2) + fv(1)) / h_d^2;
    f2v(n_dense) = (fv(n_dense) - 2*fv(n_dense-1) + fv(n_dense-2)) / h_d^2;
    for k = 2 : n_dense-1
        f2v(k) = (fv(k+1) - 2*fv(k) + fv(k-1)) / h_d^2; % centered second difference
    end
    M2(j) = max(abs(f2v));                       % keep the maximum absolute second derivative
end

% -- tolerances spanning 6 decades -- (the error targets we'll test)
tols    = [1e-2, 1e-4, 1e-6, 1e-8];              % list of tolerances to use
n_tols  = length(tols);                          % how many tolerances we have
n_funcs = 5;                                     % number of test functions

fprintf('\n============================================================\n');
fprintf('PROGRAM 4: Composite Numerical Quadrature\n');
fprintf('%d integrals, %d tolerances spanning %.0e to %.0e\n', ...
    n_funcs, n_tols, tols(1), tols(end));        % print a short summary header
fprintf('============================================================\n\n');

% ============================================================
% TASK i: A priori step size prediction for CTR
% ============================================================
% For each integral and tolerance, solve the CTR error bound for Hm:
%   |E| <= (b-a)/12 * Hm^2 * M2 <= tol  =>  Hm <= sqrt(12*tol/((b-a)*M2))
% Round m up to the next integer and recompute the actual Hm used.

fprintf('------------------------------------------------------------\n'); % print a separator line
fprintf('TASK i: A priori step size prediction (CTR)\n');              % print task title
fprintf('------------------------------------------------------------\n\n'); % print another separator and blank line

m_apriori      = zeros(n_funcs, n_tols);    % allocate array to store predicted m values
Hm_apriori     = zeros(n_funcs, n_tols);    % allocate array to store predicted Hm values
bound_apriori  = zeros(n_funcs, n_tols);    % allocate array to store the predicted error bounds

for j = 1 : n_funcs                          % loop over each test function
    a_j = limits(j, 1);                      % left endpoint of the j-th interval
    b_j = limits(j, 2);                      % right endpoint of the j-th interval
    L_j = b_j - a_j;                         % length of the integration interval

    fprintf('%s  on [%.4g, %.4g],  M_2 = %.4e\n', fnames{j}, a_j, b_j, M2(j)); % print info about this integral
    fprintf('  %-12s  %-14s  %-10s  %-14s\n', 'tolerance', 'Hm (predicted)', 'm', 'error bound'); % print column headers
    fprintf('  %s\n', repmat('-', 1, 54));    % print a line of dashes under the headers

    for it = 1 : n_tols                      % loop over each tolerance
        tau  = tols(it);                    % current tolerance value
        Hm_p = sqrt(12 * tau / (L_j * M2(j))); % compute the maximum allowed Hm from the bound
        m_p  = ceil(L_j / Hm_p);            % choose the number of subintervals by rounding up
        Hm_p = L_j / m_p;                   % recompute the actual Hm that results from integer m
        bnd  = (L_j / 12) * Hm_p^2 * M2(j); % compute the resulting error bound using the actual Hm

        m_apriori(j, it)     = m_p;         % store the predicted number of subintervals
        Hm_apriori(j, it)    = Hm_p;        % store the actual step length used
        bound_apriori(j, it) = bnd;         % store the predicted error bound

        fprintf('  %-12.2e  %-14.6e  %-10d  %-14.6e\n', tau, Hm_p, m_p, bnd); % print results for this tolerance
    end
    fprintf('\n');                           % blank line between functions
end
% ============================================================
% TASK ii: CTR with a priori Hm -- predicted vs true error
% ============================================================
% The key question here is: how tight IS the a priori bound?
% For each (integrand, tolerance) pair, compute the ratio
%   ratio = true_error / predicted_bound
% A ratio near 1 means the bound is tight. A ratio much less than 1
% means M2 bounds a worst case that the integrand rarely achieves.

fprintf('------------------------------------------------------------\n'); % print a separator line for task ii
fprintf('TASK ii: CTR with a priori Hm -- tightness of the bound\n'); % print task title
fprintf('------------------------------------------------------------\n\n'); % print separator and blank line

trap_apri_evals = zeros(n_funcs, n_tols); % initialize array to store number of function evals for a priori CTR
trap_apri_err   = zeros(n_funcs, n_tols); % initialize array to store true errors for a priori CTR

for j = 1 : n_funcs                          % loop over each test function
    a_j = limits(j, 1);                      % left endpoint of current interval
    b_j = limits(j, 2);                      % right endpoint of current interval

    fprintf('%s\n', fnames{j});              % print the name of the current function
    fprintf('  %-12s  %-10s  %-14s  %-14s  %-14s  %-10s\n', ... % print column headers
            'tolerance', 'm', 'approx Q', 'true err', 'pred bound', 'ratio');
    fprintf('  %s\n', repmat('-', 1, 76));    % print a line of dashes under the headers

    for it = 1 : n_tols                      % loop over each tolerance value
        tau = tols(it);                      % current tolerance being considered
        m_p = m_apriori(j, it);              % predicted number of subintervals for this tolerance

        [Q, nev] = comp_trap(funcs{j}, a_j, b_j, m_p); % run composite trapezoid with m_p subintervals, get approx and eval count
        true_err = abs(Q - exact(j));        % compute absolute error compared to known exact integral
        ratio    = true_err / bound_apriori(j, it); % ratio of actual error to predicted bound

        trap_apri_evals(j, it) = nev;        % store number of function evaluations used
        trap_apri_err(j, it)   = true_err;   % store the true error computed

        fprintf('  %-12.2e  %-10d  %-14.6e  %-14.6e  %-14.6e  %-10.4f\n', ... % print results for this tolerance
                tau, m_p, Q, true_err, bound_apriori(j, it), ratio);
    end
    fprintf('\n');                           % blank line between function outputs
end

% ---- Figure 1: predicted bound vs true error, log-log vs m ----
% one panel per integrand. two curves per panel: the a priori bound
% curve (from the formula) and the true error curve (from actually
% running CTR). the vertical gap between them IS the conservativeness.
% f5 should show a huge gap (M2=2000 is a worst-case corner estimate);
% f1 should show the curves nearly overlapping (f'' = e^x is smooth
% and the worst case IS the whole interval).
figure(1); clf;                              % create or clear figure 1 for plotting
sgtitle('Figure 1: CTR A Priori Bound vs True Error -- Tightness of the Bound', 'FontSize', 13); % add a title for the whole figure

for j = 1 : n_funcs                          % loop to create a subplot for each integrand
    subplot(2, 3, j);                        % place subplot in a 2x3 grid
                                             % we need a dense sweep in m to draw both curves smoothly
    a_j = limits(j, 1);                      % left endpoint for current function
    b_j = limits(j, 2);                      % right endpoint for current function
    L_j = b_j - a_j;                         % length of the integration interval

    m_sweep    = unique(round(logspace(0, 5, 30))); % choose a range of m values spaced logarithmically, make them integers
    bound_swp  = zeros(size(m_sweep));       % allocate array to hold predicted bound values for sweep
    err_swp    = zeros(size(m_sweep));       % allocate array to hold true error values for sweep
    for im = 1 : length(m_sweep)             % loop over each m in the sweep
        Hm_s         = L_j / m_sweep(im);    % compute step length corresponding to this m
        bound_swp(im) = (L_j / 12) * Hm_s^2 * M2(j); % compute predicted bound for this Hm
        [Qs, ~]      = comp_trap(funcs{j}, a_j, b_j, m_sweep(im)); % compute trapezoid approx for this m
        err_swp(im)  = max(abs(Qs - exact(j)), 1e-18); % record the true error, avoid zero for log plot
    end

    loglog(m_sweep, bound_swp, '-',  'Color', clr_bound, 'LineWidth', 2.0, ... % plot predicted bound vs m on log-log axes
           'DisplayName', 'Predicted bound');
    hold on; grid on; box on;               % keep plotting on same axes, enable grid and box
    loglog(m_sweep, err_swp,   '-',  'Color', clr_true,  'LineWidth', 2.0, ... % plot true error vs m on same axes
           'DisplayName', 'True error');
    loglog(m_apriori(j,:), trap_apri_err(j,:), 'o', 'Color', clr_apri, ... % mark the chosen m values and their true errors
           'MarkerSize', 8, 'MarkerFaceColor', clr_apri, ...
           'DisplayName', 'Chosen m(\tau)');

    xlabel('m (subintervals)'); ylabel('error'); % label axes
    title(sprintf('%s  (M_2 = %.2e)', fshort{j}, M2(j))); % title for this subplot showing short name and M2
    legend('Location', 'southwest', 'FontSize', 7); % add a legend in the lower-left of the subplot
    hold off;                               % release the current axes for the next subplot
end

% ============================================================
% TASK iii: Adaptive CTR -- compare accuracy per eval vs a priori
% ============================================================
% Task statement: "how did the accuracy per function evaluation ratio
% behave for the two sets of experiments."
%
% The key metric is the SAVINGS RATIO:
%   savings = (a priori evals) / (adaptive evals)
% Values > 1 mean adaptive needs fewer evals to hit the same tolerance.
% Higher savings correlate with looser a priori bounds (big M2, small
% average |f''|).

fprintf('------------------------------------------------------------\n'); % print separator line to console
fprintf('TASK iii: Adaptive CTR (alpha = 1/2) vs A Priori CTR\n'); % print task heading
fprintf('------------------------------------------------------------\n\n'); % print another separator and blank line

trap_adap_evals = zeros(n_funcs, n_tols); % preallocate array for adaptive trapezoid eval counts
trap_adap_err   = zeros(n_funcs, n_tols); % preallocate array for adaptive trapezoid true errors

for j = 1 : n_funcs % loop over each integrand
    a_j = limits(j, 1); % left endpoint for current integrand
    b_j = limits(j, 2); % right endpoint for current integrand

    fprintf('%s\n', fnames{j}); % print full function name
    fprintf('  %-12s  %-10s  %-10s  %-12s  %-12s  %-14s  %-14s\n', ... % print header row for the table
            'tol', 'apri m', 'adap M', 'apri evals', 'adap evals', ...
            'apri err/eval', 'adap err/eval');
    fprintf('  %s\n', repmat('-', 1, 92)); % print a line of dashes under the header

    for it = 1 : n_tols % loop over tolerances
        tau = tols(it); % current tolerance value
        [Q_ad, M_f, nev_ad] = comp_trap_adaptive(funcs{j}, a_j, b_j, tau); % run adaptive trapezoid to meet tau
        true_err_ad = abs(Q_ad - exact(j)); % compute the true absolute error for adaptive result

        trap_adap_evals(j, it) = nev_ad; % store number of function evaluations used by adaptive method
        trap_adap_err(j, it)   = true_err_ad; % store adaptive method's true error

        apri_eff = trap_apri_err(j, it) / trap_apri_evals(j, it); % compute a priori error per eval (error divided by evals)
        adap_eff = true_err_ad / nev_ad; % compute adaptive error per eval (error divided by evals)

        fprintf('  %-12.2e  %-10d  %-10d  %-12d  %-12d  %-14.4e  %-14.4e\n', ... % print a row of results for this tolerance
                tau, m_apriori(j, it), M_f, ...
                trap_apri_evals(j, it), nev_ad, apri_eff, adap_eff);
    end
    fprintf('\n'); % blank line between functions
end

% ---- Figure 2: adaptive savings ratio vs tolerance ----
% for each integrand, plot (apri_evals / adap_evals) vs tolerance.
% bars above 1 mean adaptive wins. grouped by integrand so we can
% see which ones benefit most from adaptive sizing.
% expected pattern:
%   f5 (M2 huge, f'' spikes near x=0.1): big savings
%   f4 (oscillatory, f'' nontrivial): moderate savings
%   f1 (f''=e^x smooth, M2 ~= actual): ~1
figure(2); clf; % create or clear figure 2 for plotting
sgtitle('Figure 2: Adaptive Savings -- (A Priori evals) / (Adaptive evals)', 'FontSize', 13); % add overall title for figure

savings_ratio = trap_apri_evals ./ trap_adap_evals; % compute savings ratio for each function and tolerance
bar_handle = bar(savings_ratio', 'grouped'); % create grouped bar chart with tolerances on x-axis
for j = 1 : n_funcs % loop to assign legend names for each integrand
    bar_handle(j).DisplayName = fshort{j}; % set display name for this bar group to short function name
end
hold on; grid on; box on; % keep plot for further decorations and enable grid and box
yline(1, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off'); % draw a horizontal line at 1 to show break-even
set(gca, 'XTickLabel', arrayfun(@(t) sprintf('%.0e', t), tols, 'UniformOutput', false)); % label x-ticks with tolerance values
xlabel('tolerance \tau'); ylabel('savings ratio (apri / adap)'); % label axes
title('Values > 1 mean adaptive needs fewer evals to hit the same tolerance'); % add descriptive title
legend('Location', 'northwest'); % put legend in the northwest corner
hold off; % release the plot
% ============================================================
% TASK iv: CGL vs CTR -- the O(H^2) vs O(H^4) story
% ============================================================
% Two sub-experiments:
%   A) convergence sweep: for each method, sweep M and track error.
%      confirms observed O(H^2) for both CTR variants and O(H^4) for CGL.
%   B) cost to match accuracy: for each CTR a priori tolerance, find the
%      smallest M such that CGL achieves the SAME or better true error.
%      report the ratio of CTR evals to CGL evals.

fprintf('------------------------------------------------------------\n'); % print a separator line to mark the task start
fprintf('TASK iv: CGL vs CTR -- convergence rates and cost comparison\n'); % print task title to console
fprintf('------------------------------------------------------------\n\n'); % print another separator line and a blank line

% -- sub-experiment A: convergence sweep at a common set of M values --
M_sweep = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]; % list of M values to test (number of subintervals)
nM      = length(M_sweep); % number of M values in the sweep

ctr_err_sweep = zeros(n_funcs, nM); % preallocate array to hold CTR errors for each function and M
cgl_err_sweep = zeros(n_funcs, nM); % preallocate array to hold CGL errors for each function and M

for j = 1 : n_funcs % loop over each integrand
    a_j = limits(j, 1); % left endpoint for current integrand
    b_j = limits(j, 2); % right endpoint for current integrand
    for im = 1 : nM % loop over each M in the sweep
        Mi = M_sweep(im); % current number of subintervals
        [Qt, ~] = comp_trap(funcs{j},   a_j, b_j, Mi); % compute trapezoid rule result with Mi subintervals
        [Qg, ~] = comp_gauss2(funcs{j}, a_j, b_j, Mi); % compute CGL (Gauss-Lobatto-like) result with Mi
        ctr_err_sweep(j, im) = max(abs(Qt - exact(j)), 1e-18); % store CTR error, avoid zeros by flooring
        cgl_err_sweep(j, im) = max(abs(Qg - exact(j)), 1e-18); % store CGL error, avoid zeros by flooring
    end
end

% -- sub-experiment B: CGL cost to match CTR a priori accuracy --
cgl_evals = zeros(n_funcs, n_tols); % preallocate array for CGL function evaluation counts
cgl_err   = zeros(n_funcs, n_tols); % preallocate array for CGL true errors

for j = 1 : n_funcs % loop over each integrand
    a_j = limits(j, 1); % left endpoint for current integrand
    b_j = limits(j, 2); % right endpoint for current integrand

    fprintf('%s\n', fnames{j}); % print the full name of the current function
    fprintf('  %-12s  %-14s  %-6s  %-12s  %-12s  %-14s\n', ...
            'tol', 'CTR apri err', 'CGL M', 'CGL evals', 'CTR apri evals', 'CGL true err'); % print table header
    fprintf('  %s\n', repmat('-', 1, 78)); % print a separator line under the header

    for it = 1 : n_tols % loop over each tolerance value
        target_err = trap_apri_err(j, it); % set target error equal to CTR a priori error for this tol
        M_gl = 1; % start with the smallest CGL M to try to meet the target error
        while true % loop until CGL true error <= target or an upper limit is reached
            [Qg, nev_g] = comp_gauss2(funcs{j}, a_j, b_j, M_gl); % compute CGL result and eval count for current M_gl
            terr_g      = abs(Qg - exact(j)); % compute the true absolute error for CGL
            if terr_g <= target_err || M_gl > 2^18 % stop if error met or M_gl is ridiculously large
                break
            end
            M_gl = M_gl + 1; % otherwise increase M_gl and try again
        end
        cgl_evals(j, it) = nev_g; % store number of function evaluations CGL used to meet target (or last tried)
        cgl_err(j, it)   = terr_g; % store the resulting CGL true error

        fprintf('  %-12.2e  %-14.4e  %-6d  %-12d  %-12d  %-14.4e\n', ...
                tols(it), target_err, M_gl, nev_g, trap_apri_evals(j, it), terr_g); % print a row of results
    end
    fprintf('\n'); % blank line between functions for readability
end

% ---- Figure 3: convergence rates with reference slopes ----
% log-log error vs Hm for both CTR and CGL across all integrands.
% dashed reference lines of slope 2 and 4. CTR lines should track
% slope 2; CGL lines should track slope 4. One panel per integrand.
figure(3); clf; % create or clear figure 3
sgtitle('Figure 3: Convergence Rates -- O(H^2) for CTR vs O(H^4) for CGL', 'FontSize', 13); % overall figure title

for j = 1 : n_funcs % loop to make one subplot per integrand
    subplot(2, 3, j); % position subplot in a 2-by-3 grid

    a_j = limits(j, 1); % left endpoint for current integrand
    b_j = limits(j, 2); % right endpoint for current integrand
    L_j = b_j - a_j; % interval length for current integrand
    H_sweep = L_j ./ M_sweep; % compute mesh sizes corresponding to M_sweep

    loglog(H_sweep, ctr_err_sweep(j,:), 'o-', 'Color', clr_apri, ...
           'LineWidth', 1.8, 'MarkerSize', 6, 'DisplayName', 'CTR'); % plot CTR errors vs H on log-log scale
    hold on; grid on; box on; % keep current axes and enable grid and box
    loglog(H_sweep, cgl_err_sweep(j,:), '^-', 'Color', clr_cgl, ...
           'LineWidth', 1.8, 'MarkerSize', 6, 'DisplayName', 'CGL'); % plot CGL errors vs H on log-log scale

    % reference slopes anchored to the coarsest Hm
    H_ref = logspace(log10(min(H_sweep)), log10(max(H_sweep)), 50); % make a smooth H grid for reference lines
    ref2  = ctr_err_sweep(j, 1) * (H_ref / H_sweep(1)).^2; % build an O(H^2) reference curve anchored at the first point
    ref4  = cgl_err_sweep(j, 1) * (H_ref / H_sweep(1)).^4; % build an O(H^4) reference curve anchored at the first point
    loglog(H_ref, ref2, '--', 'Color', clr_ref2, 'LineWidth', 1.2, ...
           'DisplayName', 'O(H^2) ref'); % plot the O(H^2) reference line
    loglog(H_ref, ref4, ':',  'Color', clr_ref4, 'LineWidth', 1.4, ...
           'DisplayName', 'O(H^4) ref'); % plot the O(H^4) reference line

    xlabel('H_m'); ylabel('|Q - exact|'); % label axes simply
    title(fshort{j}); % set subplot title to the short function name
    legend('Location', 'southeast', 'FontSize', 7); % add a small legend in the lower right
    hold off; % release the axes for the next subplot
end

% ---- Figure 4: error vs total function evaluations, all 3 methods ----
% this is THE efficiency figure. x-axis is total cost, y-axis is
% the error actually achieved. lower-left is better.
% CGL should dominate for tight tolerances; CTR adaptive should match
% or beat CTR apriori depending on how loose M2 is for that integrand.
figure(4); clf; % create or clear figure 4
sgtitle('Figure 4: Error vs Total Function Evaluations', 'FontSize', 13); % overall title for figure 4

for j = 1 : n_funcs % loop over each integrand to make subplots
    subplot(2, 3, j); % position subplot in a 2-by-3 layout

    loglog(trap_apri_evals(j,:), trap_apri_err(j,:), 'o-', 'Color', clr_apri, ...
           'LineWidth', 1.8, 'MarkerSize', 7, 'DisplayName', 'CTR a priori'); % plot CTR a priori error vs evals
    hold on; grid on; box on; % keep plot and enable grid and box
    loglog(trap_adap_evals(j,:), trap_adap_err(j,:), 's-', 'Color', clr_adap, ...
           'LineWidth', 1.8, 'MarkerSize', 7, 'DisplayName', 'CTR adaptive'); % plot CTR adaptive error vs evals
    loglog(cgl_evals(j,:),       cgl_err(j,:),       '^-', 'Color', clr_cgl, ...
           'LineWidth', 1.8, 'MarkerSize', 7, 'DisplayName', 'CGL'); % plot CGL error vs evals

    xlabel('total function evaluations'); ylabel('true error'); % simple axis labels
    title(fshort{j}); % set subplot title to short function name
    legend('Location', 'southwest', 'FontSize', 7); % put legend in lower-left
    hold off; % release the subplot
end

% ---- Figure 5: CGL dominance ratio vs tolerance ----
% the crossover question: for each integrand, plot the ratio
%   (CTR adaptive evals) / (CGL evals)
% across tolerances. if it's always > 1, CGL wins everywhere.
% if < 1 at loose tol and > 1 at tight tol, you see the crossover where
% CGL's O(H^4) advantage overtakes its 2x per-subinterval startup cost.
figure(5); clf; hold on; grid on; box on;                                 % make figure 5, clear it, keep plots, show grid and box
sgtitle('Figure 5: CGL Dominance -- (CTR Adaptive evals) / (CGL evals)', 'FontSize', 13); % add a big title for the whole figure

cgl_dom = trap_adap_evals ./ cgl_evals;                                     % compute ratio of CTR adaptive evals to CGL evals for each case
markers = {'o', 's', '^', 'd', 'v'};                                         % pick a set of marker symbols to cycle through
clrs_f  = lines(n_funcs);                                                    % get a distinct color for each function
for j = 1 : n_funcs                                                          % loop over each integrand/function
    semilogx(tols, cgl_dom(j, :), ['-' markers{j}], 'Color', clrs_f(j, :), ... % plot the ratio vs tolerance on a log x-axis using a unique marker and color
             'LineWidth', 1.8, 'MarkerSize', 8, 'DisplayName', fshort{j});     % set line thickness, marker size and legend entry
end                                                                           % finished plotting all functions
yline(1, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');                % draw a dashed horizontal line at 1 to show parity
set(gca, 'XDir', 'reverse');   % tighter tolerances move to the right            % reverse x-axis so tighter tolerances appear to the right
xlabel('tolerance \tau (tighter \rightarrow)'); ylabel('eval ratio (CTR adap / CGL)'); % label the axes plainly
title('Values > 1 mean CGL needs fewer evaluations; look for the crossover as \tau tightens'); % add a descriptive title for the plot
legend('Location', 'northwest');                                               % put the legend in the upper-left
hold off;                                                                      % release the current figure hold

fprintf('============================================================\n'); % print a separator line to the command window
fprintf('DONE. Produced Tables 1-4 and Figures 1-5.\n');                       % print a completion message
fprintf('============================================================\n'); % print another separator line