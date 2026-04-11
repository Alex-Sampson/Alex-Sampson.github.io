% task3_relaxation_driver.m
%
% ============================================================
% TASK 3: Scalar nonlinear relaxation (Study Set 6, Problem 6.6)
% ============================================================
%
% Iteration:
%   x_{k+1} = phi(x_k),  where
%   phi(x) = x - (x^2 - 3)/(x^2 + 2x - 3).
%
% Fixed points (solutions of x = phi(x)):
%   x =  +sqrt(3),  x = -sqrt(3)
%
% Singularities:
%   x = 1,  x = -3
%
% Figures (kept uncluttered on purpose):
%   Fig 1: Graph diagram (y=x and y=phi(x)) on a focused window + cobwebs
%   Fig 2: Residual histories for a small set of representative runs
%   Fig 3: Iterate histories (zoomed to show convergence to ±sqrt(3))
%   Fig 4: Divergent/failed iterate history (if present)
%
% Requires:
%   fixed_point.m
%
% ------------------------------------------------------------

clear; clc; close all;  % Clear the workspace, command window, and close all figures

% ---------------- user controls ----------------
tol   = 1e-12;  % Set the tolerance for convergence
maxit = 60;     % Set the maximum number of iterations allowed

den_tol = 1e-14;   % Set a threshold for the denominator to avoid singularities

alpha_plus  = sqrt(3);  % Define the positive fixed point
alpha_minus = -sqrt(3); % Define the negative fixed point

singA = -3;  % Define the first singularity
singB =  1;  % Define the second singularity

% Safe phi handle: returns NaN if too close to singularities
phi = @(x) phi_safe(x, den_tol);  % Create an anonymous function for phi that checks for singularities

fprintf('\n============================================================\n');  % Print a header for the task
fprintf('TASK 3: Scalar nonlinear relaxation (Problem 6.6)\n');  % Print the task description
fprintf('phi(x) = x - (x^2-3)/(x^2+2x-3)\n');  % Print the function phi
fprintf('fixed points: +sqrt(3)=%.16e,  -sqrt(3)=%.16e\n', alpha_plus, alpha_minus);  % Print the fixed points
fprintf('singularities: x=%.16g and x=%.16g\n', singB, singA);  % Print the singularities
fprintf('tol=%.1e, maxit=%d\n', tol, maxit);  % Print the tolerance and maximum iterations
fprintf('============================================================\n\n');  % Print a footer for the task

% ---------------- initial guesses ----------------
% Keep the plotting set SMALL so the figures are readable.
% We'll still print a summary for all runs.
cases = [ ...  % Create an array of structures for different initial guesses
    struct('x0', alpha_plus + 0.10, 'tag', 'right of +sqrt(3)', 'plot', true), ...  % Case near the positive fixed point
    struct('x0', alpha_minus - 0.10,'tag', 'left of -sqrt(3)',  'plot', true), ...  % Case near the negative fixed point
    struct('x0', 0.50,             'tag', 'inside (0, +sqrt(3))','plot', true), ...  % Case inside the interval
    struct('x0', 1.20,             'tag', 'near singularity x=1','plot', true), ...  % Case near a singularity
    struct('x0', -3.20,            'tag', 'left of x=-3',        'plot', false) ...  % Case left of the singularity
];

runs = struct();  % Initialize an empty structure to hold the results of each run

for j = 1:length(cases)  % Loop through each case
    x0 = cases(j).x0;  % Get the initial guess for this case

    [x_end, xhist, rhist, dxhist, k] = fixed_point(phi, x0, tol, maxit);  % Call the fixed point function to get results

    status = classify_run(xhist, rhist, tol, alpha_plus, alpha_minus);  % Classify the run based on the results

    runs(j).x0     = x0;  % Store the initial guess in the results structure
    runs(j).tag    = cases(j).tag;  % Store the tag for this case
    runs(j).plot   = cases(j).plot;  % Store whether to plot this case
    runs(j).x_end  = x_end;  % Store the final result of the run
    runs(j).k      = k;  % Store the number of iterations taken
    runs(j).xhist  = xhist;  % Store the history of x values
    runs(j).rhist  = rhist;  % Store the history of residuals
    runs(j).dxhist = dxhist;  % Store the history of differences
    runs(j).status = status;  % Store the status of the run

    % Short iterate print (truncate)
    fprintf('--------------------------------------------\n');  % Print a separator
    fprintf('Run %d: x0 = %+0.16e   (%s)\n', j, x0, cases(j).tag);  % Print the run number and initial guess
    fprintf('--------------------------------------------\n');  % Print another separator
    print_iter_table(xhist, rhist, [alpha_minus, alpha_plus]);  % Print the iteration table for this run

    [dist, which_root] = nearest_root(xhist(end), alpha_plus, alpha_minus);  % Find the nearest root to the final x value
    fprintf('  summary: status=%s, k=%d, x_end=%+.16e, dist_to_%s=%.3e, |r_end|=%.3e\n\n', ...  % Print a summary of the run
        status, k, xhist(end), which_root, dist, abs(rhist(end)));  % Include status, iterations, final x, distance to root, and final residual
end

% ---------------- summary table ----------------
fprintf('\n============================================================\n');  % Print a header for the summary table
fprintf('SUMMARY TABLE\n');  % Print the title of the summary table
fprintf('============================================================\n');  % Print a separator
fprintf('%-4s %-13s %-23s %-5s %-13s %-12s %-18s\n', 'run', 'x0', 'tag', 'k', 'x_end', '|r_end|', 'status');  % Print the table headers
fprintf('%s\n', repmat('-',1,90));  % Print a line separator
for j = 1:length(runs)  % Loop through each run to print the results
    fprintf('%-4d %+0.6f   %-23s %-5d %+0.6f   %-12.3e %-18s\n', ...  % Print the details of each run
        j, runs(j).x0, runs(j).tag, runs(j).k, runs(j).xhist(end), abs(runs(j).rhist(end)), runs(j).status);  % Include run number, initial guess, tag, iterations, final x, final residual, and status
end
fprintf('\n');  % Print a newline for spacing

% Choose which runs to plot (only a few)
plot_idx = find([runs.plot]);  % Find the indices of runs that should be plotted

% ============================================================
% FIGURE 1: Graph diagram + a couple cobwebs (focused window)
% ============================================================
fig1_x = [-3.5, 2.6];  % Set the x-axis limits for the figure
fig1_y = [-3.5, 3.5];  % Set the y-axis limits for the figure

figure(1); clf; hold on; grid on;  % Create a new figure and set up the grid

% y=x
xx = linspace(fig1_x(1), fig1_x(2), 1200);  % Generate x values for the line y=x
plot(xx, xx, 'k-', 'LineWidth', 1.5, 'DisplayName', 'y = x');  % Plot the line y=x

% y=phi(x) as a single polyline with NaN breaks near singularities
[xp, yp] = sample_phi_with_breaks(fig1_x(1), fig1_x(2), 4000, den_tol);  % Sample phi(x) while avoiding singularities
plot(xp, yp, '-', 'LineWidth', 1.5, 'DisplayName', 'y = \phi(x)');  % Plot the sampled phi(x)

% fixed points
plot(alpha_minus, alpha_minus, 'ko', 'MarkerFaceColor','k', 'DisplayName', 'fixed points');  % Plot the negative fixed point
plot(alpha_plus,  alpha_plus,  'ko', 'MarkerFaceColor','k', 'HandleVisibility','off');  % Plot the positive fixed point

% singularities
xline(singA, 'k--', 'LineWidth', 1.0, 'DisplayName', 'singularities');  % Draw a vertical line for the first singularity
xline(singB, 'k--', 'LineWidth', 1.0, 'HandleVisibility','off');  % Draw a vertical line for the second singularity

xlabel('x'); ylabel('y');  % Label the axes
title('Figure 1: Graph diagram for Problem 6.6 (focused window)');  % Set the title for the figure
axis([fig1_x(1) fig1_x(2) fig1_y(1) fig1_y(2)]);  % Set the axis limits

% Cobwebs: draw at most 20 steps, and stop if it leaves window
% pick two convergent (if present) and one problematic (0.5 or near singular)
idx_plus  = pick_first_with_status(runs, 'conv_to_+sqrt3');  % Find the first run that converged to the positive fixed point
idx_minus = pick_first_with_status(runs, 'conv_to_-sqrt3');  % Find the first run that converged to the negative fixed point
idx_prob  = pick_first_with_status(runs, 'fail');  % Find the first run that failed
if isnan(idx_prob)  % If no failed run was found
    % if no NaN fail, use the run labeled inside (0, +sqrt(3))
    idx_prob = pick_first_with_tag_contains(runs, 'inside');  % Use a run that is inside the interval if no failures were found
end

if ~isnan(idx_plus)  % If a convergent run to the positive fixed point was found
    add_cobweb(runs(idx_plus).xhist, phi, 20, fig1_x, fig1_y);  % Add cobwebs for this run
end
if ~isnan(idx_minus)  % If a convergent run to the negative fixed point was found
    add_cobweb(runs(idx_minus).xhist, phi, 20, fig1_x, fig1_y);  % Add cobwebs for this run
end
if ~isnan(idx_prob)  % If a problematic run was found
    add_cobweb(runs(idx_prob).xhist, phi, 12, fig1_x, fig1_y);  % Add cobwebs for this run
end

legend('Location','southeast');  % Place the legend in the southeast corner
hold off;  % Release the hold on the current figure

% ============================================================
% FIGURE 2: Error histories to the attracting fixed point
% ============================================================
figure(2); clf; hold on; grid on;  % Create a new figure for error histories, clear previous plots, and set up the grid

max_k_plot = 25;  % Set the maximum number of iterations to plot

for t = 1:length(plot_idx)  % Loop through each index in plot_idx
    j = plot_idx(t);  % Get the current index from plot_idx
    xh = runs(j).xhist;  % Extract the x history for the current run
    kk = 0:length(xh)-1;  % Create an array of iteration indices

    % Determine which fixed point the current run is approaching based on the final value
    if abs(xh(end) - alpha_plus) <= abs(xh(end) - alpha_minus)  % Check if it's closer to the positive fixed point
        alpha = alpha_plus;  % Set alpha to the positive fixed point
        lab = sprintf('x0=%+.2f \\rightarrow +\\sqrt{3}', runs(j).x0);  % Create a label for the plot
    else  % Otherwise, it's closer to the negative fixed point
        alpha = alpha_minus;  % Set alpha to the negative fixed point
        lab = sprintf('x0=%+.2f \\rightarrow -\\sqrt{3}', runs(j).x0);  % Create a label for the plot
    end

    err = abs(xh - alpha);  % Calculate the error as the absolute difference from the fixed point
    err(err==0) = eps;  % Replace any zeros in the error with a small value to avoid issues on a log scale

    keep = kk <= max_k_plot;  % Create a logical array to keep only the indices within the maximum plot limit
    kk = kk(keep);  % Filter the iteration indices to keep only those within the limit
    err = err(keep);  % Filter the error values to keep only those within the limit

    semilogy(kk, err, 'LineWidth', 1.5, 'DisplayName', lab);  % Plot the error on a semi-logarithmic scale
end

set(gca,'YScale','log');  % Set the y-axis to a logarithmic scale
xlabel('k');  % Label the x-axis
ylabel('|x_k - \alpha|');  % Label the y-axis with the error
title(sprintf('Figure 2: Error histories (k \\le %d)', max_k_plot));  % Set the title for the figure
legend('Location','eastoutside');  % Position the legend outside the plot area on the east side
xlim([0, max_k_plot]);  % Set the x-axis limits
grid on; grid minor;  % Enable grid lines for better readability
hold off;  % Release the hold on the current figure

% ============================================================
% FIGURE 3: Iterate histories (zoomed to show convergence)
% ============================================================
figure(3); clf; hold on; grid on;  % Create a new figure for iterate histories, clear previous plots, and set up the grid

for t = 1:length(plot_idx)  % Loop through each index in plot_idx
    j = plot_idx(t);  % Get the current index from plot_idx
    xh = runs(j).xhist;  % Extract the x history for the current run
    name = legend_name(runs(j), alpha_plus, alpha_minus);  % Get the legend name for the current run
    plot(0:length(xh)-1, xh, '-o', 'LineWidth', 1.0, 'MarkerSize', 3, 'DisplayName', name);  % Plot the iterate history with markers
end

yline(alpha_minus, 'k--', 'LineWidth', 1.0, 'DisplayName', '-\sqrt{3}');  % Draw a horizontal line for the negative fixed point
yline(alpha_plus,  'k--', 'LineWidth', 1.0, 'DisplayName', '+\sqrt{3}');  % Draw a horizontal line for the positive fixed point
yline(singA, 'k:', 'LineWidth', 1.0, 'DisplayName', 'x=-3');  % Draw a dashed line for the first singularity
yline(singB, 'k:', 'LineWidth', 1.0, 'DisplayName', 'x=1');  % Draw a dashed line for the second singularity

xlabel('k');  % Label the x-axis
ylabel('x_k');  % Label the y-axis with the iterate values
title('Figure 3: Iterates vs k (zoomed)');  % Set the title for the figure
legend('Location','eastoutside');  % Position the legend outside the plot area on the east side

xlim([0, 15]);  % Set the x-axis limits to zoom in on the convergence
ylim([-4, 4]);  % Set the y-axis limits to focus on the relevant range
hold off;  % Release the hold on the current figure

% ============================================================
% FIGURE 4: Nonconvergent run diagnostics (drift)
% ============================================================
idx_bad = pick_first_with_status(runs, 'maxit_or_stalled');  % Find the first nonconvergent run
if isnan(idx_bad), idx_bad = pick_first_with_status(runs, 'fail'); end  % If none found, look for a failed run

if ~isnan(idx_bad)  % If a problematic run was found
    xh  = runs(idx_bad).xhist;  % Extract the x history for the problematic run
    dxh = diff(xh);  % Calculate the differences between consecutive iterates
    rh  = abs(runs(idx_bad).rhist);  % Get the absolute values of the residuals for the problematic run

    kx  = 0:length(xh)-1;  % Create an array of iteration indices
    kdx = 0:length(dxh)-1;  % Create an array of indices for the differences

    rh(rh==0) = eps;  % Replace any zeros in the residuals with a small value to avoid issues

    figure(4); clf;  % Create a new figure for nonconvergent run diagnostics and clear previous plots

    subplot(3,1,1);  % Create the first subplot for the iterate history
    plot(kx, xh, '-o', 'LineWidth', 1.25, 'MarkerSize', 3);  % Plot the iterate history with markers
    grid on;  % Enable grid lines for better readability
    xlabel('k'); ylabel('x_k');  % Label the axes
    title(sprintf('Figure 4: Nonconvergent run (x0=%+.2f)', runs(idx_bad).x0));  % Set the title for the subplot

    subplot(3,1,2);  % Create the second subplot for the differences
    plot(kdx, dxh, '-o', 'LineWidth', 1.25, 'MarkerSize', 3);  % Plot the differences with markers
    grid on;  % Enable grid lines for better readability
    xlabel('k'); ylabel('dx_k = x_{k+1}-x_k');  % Label the axes
    yline(-1,'k--','dx=-1 (drift ref)');  % Draw a reference line for drift

    subplot(3,1,3);  % Create the third subplot for the absolute differences
    plot(kdx, abs(dxh), '-o', 'LineWidth', 1.25, 'MarkerSize', 3);  % Plot the absolute differences with markers
    grid on;  % Enable grid lines for better readability
    xlabel('k'); ylabel('|dx_k| (this equals |r_k| here)');  % Label the axes
    yline(1,'k--','|dx|=1 ref');  % Draw a reference line for the absolute difference
end

% ========================================================================
% Local helper functions (MATLAB allows functions at end of script)
% ========================================================================

function y = phi_safe(x, den_tol)
    % Evaluate phi(x) and return NaN if too close to singularities.
    den = x.^2 + 2*x - 3;  % Calculate the denominator (x-1)(x+3)
    if ~isfinite(x) || abs(den) < den_tol  % Check if x is finite and if the denominator is too small
        y = NaN;  % If so, return NaN to avoid issues
        return;  % Exit the function
    end
    y = x - (x.^2 - 3)./den;  % Compute the value of phi(x)
    if ~isfinite(y)  % Check if the result is finite
        y = NaN;  % If not, return NaN
    end
end

function [xs, ys] = sample_phi_with_breaks(a, b, n, den_tol)
    % Sample phi(x) on [a,b] and insert NaNs near singularities so the plot
    % doesn't draw vertical spikes across poles.
    xs = linspace(a, b, n);  % Create n evenly spaced points between a and b
    ys = NaN(size(xs));  % Initialize the output array with NaNs

    for i = 1:n  % Loop through each sample point
        x = xs(i);  % Get the current sample point
        den = x.^2 + 2*x - 3;  % Calculate the denominator for the current x
        if abs(den) < den_tol  % Check if the denominator is too small
            ys(i) = NaN;  % If so, set the corresponding y value to NaN
        else
            y = x - (x.^2 - 3)./den;  % Compute the value of phi(x)
            % If y is huge, it will wreck the visible scale; break instead.
            if ~isfinite(y) || abs(y) > 50  % Check if y is finite and within a reasonable range
                ys(i) = NaN;  % If not, set the corresponding y value to NaN
            else
                ys(i) = y;  % Otherwise, store the computed value
            end
        end
    end
end

function status = classify_run(xhist, rhist, tol, alpha_plus, alpha_minus)
    x_end = xhist(end);  % Get the last value of the x history
    r_end = rhist(end);  % Get the last value of the residual history

    if ~isfinite(x_end) || ~isfinite(r_end) || isnan(x_end) || isnan(r_end)  % Check if the last values are valid
        status = 'fail';  % If not, classify the run as a failure
        return;  % Exit the function
    end

    if abs(r_end) < tol  % Check if the last residual is within the tolerance
        if abs(x_end - alpha_plus) <= abs(x_end - alpha_minus)  % Determine which fixed point it's closer to
            status = 'conv_to_+sqrt3';  % Classify as converging to the positive fixed point
        else
            status = 'conv_to_-sqrt3';  % Classify as converging to the negative fixed point
        end
        return;  % Exit the function
    end

    % If it didn't converge, but ended close to a root, mark as "near".
    if min(abs(x_end - alpha_plus), abs(x_end - alpha_minus)) < 1e-8  % Check if it's near a root
        if abs(x_end - alpha_plus) <= abs(x_end - alpha_minus)  % Determine which root it's closer to
            status = 'near_+sqrt3';  % Classify as near the positive fixed point
        else
            status = 'near_-sqrt3';  % Classify as near the negative fixed point
        end
        return;  % Exit the function
    end

    status = 'maxit_or_stalled';  % If none of the above, classify as stalled or max iterations reached
end

function print_iter_table(xhist, rhist, roots)
    % Print k, x_k, |r_k|, dist to nearest root.
    K = length(xhist)-1;  % Get the number of iterations
    max_lines = 10;  % Set a limit for how many lines to print

    fprintf('  %-4s %-22s %-12s %-18s\n', 'k', 'x_k', '|r_k|', 'dist_to_nearest_root');  % Print header
    fprintf('  %s\n', repmat('-',1,62));  % Print a separator line

    for kk = 0:K  % Loop through each iteration
        xk = xhist(kk+1);  % Get the current x value
        rk = rhist(kk+1);  % Get the current residual value
        dist = min(abs(xk - roots(1)), abs(xk - roots(2)));  % Calculate distance to the nearest root

        if kk <= max_lines-1 || kk == K  % Check if we are within the line limit or at the last line
            fprintf('  %-4d %+0.16e %-12.3e %-18.3e\n', kk, xk, abs(rk), dist);  % Print the current iteration details
        elseif kk == max_lines  % If we've reached the max lines
            fprintf('  ...\n');  % Indicate that there are more lines not shown
        end
    end
end

function [dist, which_root] = nearest_root(x, alpha_plus, alpha_minus)
    dplus  = abs(x - alpha_plus);  % Calculate distance to the positive fixed point
    dminus = abs(x - alpha_minus);  % Calculate distance to the negative fixed point
    if dplus <= dminus  % Check which distance is smaller
        dist = dplus;  % Set the distance to the positive fixed point
        which_root = '+sqrt3';  % Indicate that it's the positive fixed point
    else
        dist = dminus;  % Set the distance to the negative fixed point
        which_root = '-sqrt3';  % Indicate that it's the negative fixed point
    end
end

function name = legend_name(run, alpha_plus, alpha_minus)
    % Short legend strings to avoid a giant unreadable legend.
    x0 = run.x0;  % Get the initial value for the run
    st = run.status;  % Get the status of the run

    if strcmp(st,'conv_to_+sqrt3')  % Check if the run converged to the positive fixed point
        name = sprintf('x0=%+.2f \\rightarrow +\\sqrt{3}', x0);  % Create a legend entry for this case
    elseif strcmp(st,'conv_to_-sqrt3')  % Check if the run converged to the negative fixed point
        name = sprintf('x0=%+.2f \\rightarrow -\\sqrt{3}', x0);  % Create a legend entry for this case
    elseif strcmp(st,'fail')  % Check if the run failed
        name = sprintf('x0=%+.2f (fail)', x0);  % Create a legend entry indicating failure
    else
        % classify by nearest root if meaningful
        x_end = run.xhist(end);  % Get the last x value from the history
        if isfinite(x_end)  % Check if the last x value is finite
            if abs(x_end-alpha_plus) <= abs(x_end-alpha_minus)  % Determine which root it's closer to
                name = sprintf('x0=%+.2f (stalled near +\\sqrt{3}?)', x0);  % Create a legend entry for being near the positive root
            else
                name = sprintf('x0=%+.2f (stalled near -\\sqrt{3}?)', x0);  % Create a legend entry for being near the negative root
            end
        else
            name = sprintf('x0=%+.2f (stalled)', x0);  % Create a legend entry indicating stall
        end
    end
end

function add_cobweb(xhist, phi, max_steps, xwin, ywin)
    % Draw a cobweb using the stored xhist.
    % Stops early if it leaves the plotting window.

    K = min(length(xhist)-1, max_steps);  % Determine how many steps to plot, limited by max_steps

    for k = 1:K  % Loop through each step
        xk   = xhist(k);  % Get the current x value
        xkp1 = xhist(k+1);  % Get the next x value
        if ~isfinite(xk) || ~isfinite(xkp1)  % Check if the current or next x value is valid
            return;  % If not, exit the function
        end

        yk = phi(xk);  % Compute the corresponding y value using phi
        if ~isfinite(yk)  % Check if the y value is valid
            return;  % If not, exit the function
        end

        % stop if we leave the visible window
        if xk < xwin(1) || xk > xwin(2) || yk < ywin(1) || yk > ywin(2)  % Check if we're outside the plotting window
            return;  % If so, exit the function
        end

        % vertical: (xk, xk) -> (xk, phi(xk))
        plot([xk xk], [xk yk], 'b-', 'LineWidth', 1.0, 'HandleVisibility','off');  % Draw the vertical line in the cobweb

        % horizontal: (xk, phi(xk)) -> (phi(xk), phi(xk))
        plot([xk yk], [yk yk], 'b-', 'LineWidth', 1.0, 'HandleVisibility','off');  % Draw the horizontal line in the cobweb
    end
end

function idx = pick_first_with_status(runs, wanted)
    idx = NaN;  % Initialize index as NaN
    for j = 1:length(runs)  % Loop through each run
        if strcmp(runs(j).status, wanted)  % Check if the status matches the wanted status
            idx = j;  % If it does, set the index to the current run
            return;  % Exit the function
        end
    end
end

function idx = pick_first_with_tag_contains(runs, pat)
    idx = NaN;  % Initialize index as NaN
    for j = 1:length(runs)  % Loop through each run
        if contains(runs(j).tag, pat)  % Check if the tag contains the specified pattern
            idx = j;  % If it does, set the index to the current run
            return;  % Exit the function
        end
    end
end
