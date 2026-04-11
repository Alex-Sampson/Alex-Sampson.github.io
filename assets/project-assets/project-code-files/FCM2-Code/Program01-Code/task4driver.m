% Task 4: Nonlinear relaxation for systems (Problem 6.10)
%
% System:
%   xi^2 + eta^2 = 4
%   exp(xi) + eta = 1   (eta = 1 - exp(xi))
%
% Fixed-point maps:
%   G1(x) = ( log(1-eta),    -sqrt(4 - xi^2) )
%   G2(x) = ( -sqrt(4-eta^2), 1 - exp(xi)    )
%
% What this driver produces (simple & interpretable):
%   Fig 1: Geometry + both intersections
%   Fig 2: G1 trajectories (few x0)
%   Fig 3: G1 residual histories (convergence evidence)
%   Fig 4: G2 trajectories (few x0)
%   Fig 5: G2 residual histories (convergence evidence)
%   Fig 6: Grid robustness summary (counts + iteration histogram)
%
% Requires in same folder:
%   G1.m, G2.m, fixed_point_vec.m, newton.m

clear; clc; close all;

tol   = 1e-12;
maxit = 60;

% ------------------------------------------------------------
% 1) Compute the two true intersections (reference solutions)
% ------------------------------------------------------------
% Reduce to scalar equation in xi:
%   eta = 1 - exp(xi)
%   circle => xi^2 + (1-exp(xi))^2 = 4
h  = @(xi) xi.^2 + (1 - exp(xi)).^2 - 4;  % Define the function h as a function of xi, representing the equation of the circle minus 4
hp = @(xi) 2*xi + 2*(1 - exp(xi)).*(-exp(xi));  % Define the derivative of h, hp, using the chain rule

[xi_pos,~,~,~,~] = newton(h, hp,  0.8, 1e-14, 200);  % Use Newton's method to find the positive root, starting from 0.8
[xi_neg,~,~,~,~] = newton(h, hp, -1.2, 1e-14, 200);  % Use Newton's method to find the negative root, starting from -1.2

eta_pos = 1 - exp(xi_pos);  % Calculate eta for the positive root using the equation eta = 1 - exp(xi)
eta_neg = 1 - exp(xi_neg);  % Calculate eta for the negative root using the same equation

root_pos = [xi_pos; eta_pos];  % Create a vector for the positive root (xi>0, eta<0)
root_neg = [xi_neg; eta_neg];  % Create a vector for the negative root (xi<0, eta>0)

fprintf('============================================================\n');  % Print a separator line for clarity
fprintf('TASK 4: Nonlinear relaxation for systems (Problem 6.10)\n');  % Print the task title
fprintf('tol = %.1e, maxit = %d\n', tol, maxit);  % Print the tolerance and maximum iterations used
fprintf('============================================================\n\n');  % Print another separator line

fprintf('Reference intersections (solutions of the system):\n');  % Indicate that the following lines will show the reference solutions
fprintf('  root_pos (xi>0, eta<0) = (%.16e, %.16e)\n', root_pos(1), root_pos(2));  % Print the positive root with high precision
fprintf('  root_neg (xi<0, eta>0) = (%.16e, %.16e)\n\n', root_neg(1), root_neg(2));  % Print the negative root with high precision

% ------------------------------------------------------------
% 2) Figure 1: geometry
% ------------------------------------------------------------
figure(1); clf; hold on; grid on;  % Create a new figure, clear it, hold on for multiple plots, and enable the grid
plot_system_curves();  % Plot the system curves defined in the function plot_system_curves
plot(root_pos(1), root_pos(2), 'ko', 'MarkerSize', 7, 'LineWidth', 1.5, 'DisplayName','root\_pos');  % Plot the positive root as a black circle with specified size and line width
plot(root_neg(1), root_neg(2), 'ks', 'MarkerSize', 7, 'LineWidth', 1.5, 'DisplayName','root\_neg');  % Plot the negative root as a black square with specified size and line width

xlabel('\xi');  % Label the x-axis with the symbol xi
ylabel('\eta');  % Label the y-axis with the symbol eta
title('Geometry: \xi^2+\eta^2=4 and \eta=1-e^\xi (two intersections)');  % Set the title of the plot to describe the geometry being shown
legend('Location','best');  % Add a legend to the plot, placing it in the best location
xlim([-2 2]);  % Set the limits of the x-axis from -2 to 2
ylim([-2.5 2]);  % Set the limits of the y-axis from -2.5 to 2

% ------------------------------------------------------------
% 3) Representative x0 (few, but we ALSO do a grid sweep later)
% ------------------------------------------------------------
rep = { ...  % Create a cell array named 'rep' to hold different initial guess vectors for the fixed-point iterations
    [ 0.50; -0.50], ...  % First initial guess: a vector with xi = 0.50 and eta = -0.50
    [ 1.20; -1.50], ...  % Second initial guess: a vector with xi = 1.20 and eta = -1.50
    [-0.50;  0.50], ...  % Third initial guess: a vector with xi = -0.50 and eta = 0.50
    [-1.20;  1.20]  ...  % Fourth initial guess: a vector with xi = -1.20 and eta = 1.20
};

% ------------------------------------------------------------
% 4) Figures 23: G1 trajectories + residual histories
% ------------------------------------------------------------
figure(2); clf; hold on; grid on;  % Create figure 2, clear it, hold on for multiple plots, and enable the grid
plot_system_curves();  % Call the function to plot the system curves
plot(root_pos(1), root_pos(2), 'ko', 'MarkerSize', 6, 'LineWidth', 1.2, 'DisplayName','root\_pos');  % Plot the positive root as a black circle
plot(root_neg(1), root_neg(2), 'ks', 'MarkerSize', 6, 'LineWidth', 1.2, 'DisplayName','root\_neg');  % Plot the negative root as a black square
xlabel('\xi');  % Label the x-axis with the symbol xi
ylabel('\eta');  % Label the y-axis with the symbol eta
title('G1 trajectories in (\xi,\eta)-plane');  % Set the title for the plot indicating it's for G1 trajectories
xlim([-2 2]);  % Set the limits of the x-axis from -2 to 2
ylim([-2.5 2]);  % Set the limits of the y-axis from -2.5 to 2

figure(3); clf; hold on; grid on;  % Create figure 3, clear it, hold on for multiple plots, and enable the grid
title('G1 convergence evidence: r_k = ||x_k - G1(x_k)||_2');  % Set the title for the convergence evidence plot
xlabel('k');  % Label the x-axis with 'k' for iterations
ylabel('r_k');  % Label the y-axis with 'r_k' for residuals

fprintf('-------------------------------\n');  % Print a separator line for clarity
fprintf('Representative runs for G1\n');  % Indicate that the following lines will show representative runs for G1
fprintf('Convergence test used: isfinite(x_end), isfinite(G(x_end)), ||x_end-G(x_end)||_2 < tol\n');  % Describe the convergence test used
fprintf('-------------------------------\n');  % Print another separator line

for i = 1:numel(rep)  % Loop over each initial guess in the rep cell array
    x0 = rep{i};  % Get the current initial guess

    [xend, xhist, rhist, ~, k] = fixed_point_vec(@G1, x0, tol, maxit);  % Perform fixed-point iteration and get results
    [success, r_end] = final_converged(@G1, xend, tol);  % Check if the iteration converged

    % trajectory plot
    figure(2);  % Switch to figure 2 for trajectory plotting
    plot_trajectory(xhist, sprintf('x0=(%.2f,%.2f)', x0(1), x0(2)));  % Plot the trajectory of the current initial guess

    % residual history plot
    figure(3);  % Switch to figure 3 for residual history plotting
    if ~isempty(rhist) && all(isfinite(rhist))  % Check if the residual history is not empty and all values are finite
        rplot = max(rhist(:), eps);  % Get the maximum residual, ensuring it's at least epsilon to avoid log(0)
        semilogy(0:numel(rplot)-1, rplot, '-o', 'MarkerSize', 3, ...  % Plot the residuals on a logarithmic scale
            'DisplayName', sprintf('x0=(%.2f,%.2f)', x0(1), x0(2)));  % Add a label for the current initial guess
    end

    % print summary
    if success  % If the iteration was successful
        dpos = norm(xend - root_pos);  % Calculate the distance to the positive root
        dneg = norm(xend - root_neg);  % Calculate the distance to the negative root
        fprintf('x0=(%+.2f,%+.2f): CONV  k=%2d  r_end=%.2e  dist(pos)=%.2e  dist(neg)=%.2e\n', ...  % Print a summary of the results
            x0(1), x0(2), k, r_end, dpos, dneg);
    else  % If the iteration failed
        fprintf('x0=(%+.2f,%+.2f): FAIL  k=%2d  (last computed r ~= %s)\n', ...  % Print a failure message
            x0(1), x0(2), k, num2str(r_end,'%.2e'));
    end
end

figure(2); legend('Location','bestoutside');  % Add a legend to figure 2, placing it outside the best location
figure(3); legend('Location','bestoutside');  % Add a legend to figure 3, placing it outside the best location

% ------------------------------------------------------------
% 5) Figures 45: G2 trajectories + residual histories
% ------------------------------------------------------------
figure(4); clf; hold on; grid on;  % Create figure 4, clear it, hold on for multiple plots, and enable the grid
plot_system_curves();  % Plot the system curves for G2
plot(root_pos(1), root_pos(2), 'ko', 'MarkerSize', 6, 'LineWidth', 1.2, 'DisplayName','root\_pos');  % Plot the positive root as a black circle
plot(root_neg(1), root_neg(2), 'ks', 'MarkerSize', 6, 'LineWidth', 1.2, 'DisplayName','root\_neg');  % Plot the negative root as a black square
xlabel('\xi'); ylabel('\eta');  % Label the x-axis as '\xi' and the y-axis as '\eta'
title('G2 trajectories in (\xi,\eta)-plane');  % Set the title for the G2 trajectories plot
xlim([-2 2]); ylim([-2.5 2]);  % Set the limits for the x-axis and y-axis

figure(5); clf; hold on; grid on;  % Create figure 5, clear it, hold on for multiple plots, and enable the grid
title('G2 convergence evidence: r_k = ||x_k - G2(x_k)||_2');  % Set the title for the convergence evidence plot
xlabel('k'); ylabel('r_k');  % Label the x-axis with 'k' for iterations and the y-axis with 'r_k' for residuals

fprintf('\n-------------------------------\n');  % Print a separator line for clarity
fprintf('Representative runs for G2\n');  % Indicate that the following lines will show representative runs for G2
fprintf('Convergence test used: isfinite(x_end), isfinite(G(x_end)), ||x_end-G(x_end)||_2 < tol\n');  % Describe the convergence test used
fprintf('-------------------------------\n');  % Print another separator line

for i = 1:numel(rep)  % Loop over each initial guess in the rep cell array
    x0 = rep{i};  % Get the current initial guess

    [xend, xhist, rhist, ~, k] = fixed_point_vec(@G2, x0, tol, maxit);  % Perform fixed-point iteration and get results
    [success, r_end] = final_converged(@G2, xend, tol);  % Check if the iteration converged

    % trajectory plot
    figure(4);  % Switch to figure 4 for trajectory plotting
    plot_trajectory(xhist, sprintf('x0=(%.2f,%.2f)', x0(1), x0(2)));  % Plot the trajectory of the current initial guess

    % residual history plot
    figure(5);  % Switch to figure 5 for residual history plotting
    if ~isempty(rhist) && all(isfinite(rhist) )  % Check if the residual history is not empty and all values are finite
        rplot = max(rhist(:), eps);  % Get the maximum residual, ensuring it's at least epsilon to avoid log(0)
        semilogy(0:numel(rplot)-1, rplot, '-o', 'MarkerSize', 3, ...  % Plot the residuals on a logarithmic scale
            'DisplayName', sprintf('x0=(%.2f,%.2f)', x0(1), x0(2)));  % Add a label for the current initial guess
    end

    % print summary
    if success  % If the iteration was successful
        dpos = norm(xend - root_pos);  % Calculate the distance to the positive root
        dneg = norm(xend - root_neg);  % Calculate the distance to the negative root
        fprintf('x0=(%+.2f,%+.2f): CONV  k=%2d  r_end=%.2e  dist(pos)=%.2e  dist(neg)=%.2e\n', ...  % Print a summary of the results
            x0(1), x0(2), k, r_end, dpos, dneg);
    else  % If the iteration failed
        fprintf('x0=(%+.2f,%+.2f): FAIL  k=%2d  (last computed r ~= %s)\n', ...  % Print a failure message
            x0(1), x0(2), k, num2str(r_end,'%.2e'));
    end
end

figure(4); legend('Location','bestoutside');  % Add a legend to figure 4, placing it outside the best location
figure(5); legend('Location','bestoutside');  % Add a legend to figure 5, placing it outside the best location

% ------------------------------------------------------------
% 6) Robustness: grid sweep summary (NO basins)
% ------------------------------------------------------------
% We report:
%   - how many starts are "valid" (first step defined/real)
%   - how many converge
%   - how many fail
%   - among converged: which root, and iteration statistics

xi_vals  = linspace(-2, 2, 61);  % Create a vector of 61 linearly spaced values from -2 to 2 for xi
eta_vals = linspace(-2, 2, 61);  % Create a vector of 61 linearly spaced values from -2 to 2 for eta

root_tol = 1e-7;  % Set a tolerance level for determining how close we are to the roots

% --- G1 counters ---
G1_valid = 0; G1_conv = 0; G1_fail = 0;  % Initialize counters for G1: valid starts, converged, and failed
G1_to_pos = 0; G1_to_neg = 0;  % Initialize counters for how many converge to positive and negative roots
G1_k_list = [];  % Initialize a list to store iteration counts for G1

% --- G2 counters ---
G2_valid = 0; G2_conv = 0; G2_fail = 0;  % Initialize counters for G2: valid starts, converged, and failed
G2_to_pos = 0; G2_to_neg = 0;  % Initialize counters for how many converge to positive and negative roots
G2_k_list = [];  % Initialize a list to store iteration counts for G2

for xi0 = xi_vals  % Loop over each value of xi
    for eta0 = eta_vals  % Loop over each value of eta
        x0 = [xi0; eta0];  % Create a column vector for the current initial guess

        % ---- G1 ----
        y0 = G1(x0);  % Evaluate G1 at the current guess to check if the first step is valid (real domain)
        if all(isfinite(y0))  % Check if all elements of y0 are finite
            G1_valid = G1_valid + 1;  % Increment the count of valid starts for G1
            [xend, ~, ~, ~, k] = fixed_point_vec(@G1, x0, tol, maxit);  % Perform fixed-point iteration for G1
            [success, ~] = final_converged(@G1, xend, tol);  % Check if the iteration converged successfully

            if success  % If the iteration was successful
                G1_conv = G1_conv + 1;  % Increment the count of converged runs for G1
                G1_k_list(end+1,1) = k;  % Store the number of iterations taken to converge

                if norm(xend - root_pos) < root_tol  % Check if the result is close to the positive root
                    G1_to_pos = G1_to_pos + 1;  % Increment the count for convergence to the positive root
                elseif norm(xend - root_neg) < root_tol  % Check if the result is close to the negative root
                    G1_to_neg = G1_to_neg + 1;  % Increment the count for convergence to the negative root
                else  % If it converged but not close to either root
                    G1_fail = G1_fail + 1;  % Count this as a failure for interpretation
                end
            else  % If the iteration failed
                G1_fail = G1_fail + 1;  % Increment the count of failures for G1
            end
        end

        % ---- G2 ----
        y0 = G2(x0);  % Evaluate G2 at the current guess to check if the first step is valid (real domain)
        if all(isfinite(y0))  % Check if all elements of y0 are finite
            G2_valid = G2_valid + 1;  % Increment the count of valid starts for G2
            [xend, ~, ~, ~, k] = fixed_point_vec(@G2, x0, tol, maxit);  % Perform fixed-point iteration for G2
            [success, ~] = final_converged(@G2, xend, tol);  % Check if the iteration converged successfully

            if success  % If the iteration was successful
                G2_conv = G2_conv + 1;  % Increment the count of converged runs for G2
                G2_k_list(end+1,1) = k;   % Store the number of iterations taken to converge

                if norm(xend - root_pos) < root_tol  % Check if the result is close to the positive root
                    G2_to_pos = G2_to_pos + 1;  % Increment the count for convergence to the positive root
                elseif norm(xend - root_neg) < root_tol  % Check if the result is close to the negative root
                    G2_to_neg = G2_to_neg + 1;  % Increment the count for convergence to the negative root
                else  % If it converged but not close to either root
                    G2_fail = G2_fail + 1;  % Count this as a failure for interpretation
                end
            else  % If the iteration failed
                G2_fail = G2_fail + 1;  % Increment the count of failures for G2
            end
        end
    end
end

fprintf('\n============================================================\n'); % Print a separator line for clarity in the output
fprintf('Grid robustness summary (61x61 starts in [-2,2]^2)\n'); % Indicate the summary of the grid robustness test
fprintf('Valid start = first step finite (G(x0) real)\n'); % Explain what constitutes a valid start
fprintf('Converged   = final residual check at x_end: ||x_end - G(x_end)||_2 < tol\n'); % Define what it means for the iteration to have converged
fprintf('============================================================\n'); % Print another separator line

fprintf('\nG1:\n'); % Start the summary for G1
fprintf('  valid starts: %d\n', G1_valid); % Display the count of valid starts for G1
fprintf('  converged:    %d\n', G1_conv); % Show how many iterations converged for G1
fprintf('  failed:       %d\n', G1_fail); % Indicate how many iterations failed for G1
fprintf('  to root_pos (xi>0,eta<0): %d\n', G1_to_pos); % Count of converged iterations to the positive root for G1
fprintf('  to root_neg (xi<0,eta>0): %d\n', G1_to_neg); % Count of converged iterations to the negative root for G1

if ~isempty(G1_k_list) % Check if there are any iterations recorded for G1
    fprintf('  k stats (converged only): min=%d, median=%d, max=%d, mean=%.2f\n', ... % Print statistics for the number of iterations taken to converge
        min(G1_k_list), median(G1_k_list), max(G1_k_list), mean(G1_k_list)); % Calculate and display min, median, max, and mean of k
end

fprintf('\nG2:\n'); % Start the summary for G2
fprintf('  valid starts: %d\n', G2_valid); % Display the count of valid starts for G2
fprintf('  converged:    %d\n', G2_conv); % Show how many iterations converged for G2
fprintf('  failed:       %d\n', G2_fail); % Indicate how many iterations failed for G2
fprintf('  to root_pos (xi>0,eta<0): %d\n', G2_to_pos); % Count of converged iterations to the positive root for G2
fprintf('  to root_neg (xi<0,eta>0): %d\n', G2_to_neg); % Count of converged iterations to the negative root for G2

if ~isempty(G2_k_list) % Check if there are any iterations recorded for G2
    fprintf('  k stats (converged only): min=%d, median=%d, max=%d, mean=%.2f\n', ... % Print statistics for the number of iterations taken to converge
        min(G2_k_list), median(G2_k_list), max(G2_k_list), mean(G2_k_list)); % Calculate and display min, median, max, and mean of k
end

% Figure 6: simple robustness visualization (counts + iteration hist)
figure(6); clf; % Create a new figure for visualizing the results and clear any previous content

subplot(1,2,1); % Create the first subplot in a 1x2 grid
bar([G1_conv, G1_fail; G2_conv, G2_fail]); % Create a bar graph showing converged and failed counts for G1 and G2
grid on; % Turn on the grid for better readability
set(gca,'XTickLabel',{'G1','G2'}); % Set the x-axis labels to G1 and G2
legend({'converged','failed'}, 'Location','best'); % Add a legend to the bar graph
title('Grid sweep outcomes (valid starts only)'); % Set the title for the first subplot
ylabel('count'); % Label the y-axis

subplot(1,2,2); % Create the second subplot in a 1x2 grid
hold on; grid on; % Hold the current plot and turn on the grid
if ~isempty(G1_k_list) % Check if there are any iterations recorded for G1
    histogram(G1_k_list, 'DisplayName','G1 (k for converged)'); % Create a histogram for the number of iterations for G1
end
if ~isempty(G2_k_list) % Check if there are any iterations recorded for G2
    histogram(G2_k_list, 'DisplayName','G2 (k for converged)'); % Create a histogram for the number of iterations for G2
end
title('Iteration counts for converged runs'); % Set the title for the second subplot
xlabel('k'); ylabel('frequency'); % Label the x-axis and y-axis
legend('Location','best'); % Add a legend to the histogram

fprintf('\nDone. Figures 1-6 created.\n'); % Print a completion message indicating that the figures have been created

% ============================================================
% Local helpers
% ============================================================

function plot_system_curves()  % Define a function to plot the system curves
    th = linspace(0, 2*pi, 800);  % Create a vector of angles from 0 to 2*pi with 800 points
    plot(2*cos(th), 2*sin(th), 'k-', 'LineWidth', 1.25, ...  % Plot a circle with radius 2
        'DisplayName','\xi^2+\eta^2=4');  % Label the plot for the legend

    xi = linspace(-2, 2, 800);  % Create a vector of x values from -2 to 2 with 800 points
    eta = 1 - exp(xi);  % Calculate the corresponding y values using the equation eta = 1 - e^xi
    plot(xi, eta, 'k--', 'LineWidth', 1.25, ...  % Plot the curve defined by the equation
        'DisplayName','\eta=1-e^\xi');  % Label the plot for the legend
end

function plot_trajectory(xhist, labelstr)  % Define a function to plot the trajectory of points
    if isempty(xhist), return; end  % If the history is empty, just exit the function
    good = all(isfinite(xhist), 1);  % Check if all points in the history are finite
    xhist = xhist(:, good);  % Keep only the finite points
    if size(xhist,2) < 2, return; end  % If there are less than 2 points, exit the function
    plot(xhist(1,:), xhist(2,:), '-o', 'LineWidth', 1.0, 'MarkerSize', 3, ...  % Plot the trajectory with markers
        'DisplayName', labelstr);  % Use the provided label for the legend
end

function [success, r_end] = final_converged(G, xend, tol)  % Define a function to check for convergence
    % Strong convergence check at the returned endpoint:
    %   require xend finite AND G(xend) finite AND residual < tol
    if any(~isfinite(xend))  % Check if the endpoint is finite
        success = false;  % If not, set success to false
        r_end = NaN;  % Set the residual to NaN
        return;  % Exit the function
    end
    Gx = G(xend);  % Evaluate G at the endpoint
    if any(~isfinite(Gx))  % Check if the result of G is finite
        success = false;  % If not, set success to false
        r_end = NaN;  % Set the residual to NaN
        return;  % Exit the function
    end
    r_end = norm(xend - Gx);  % Calculate the residual as the norm of the difference
    success = (r_end < tol);  % Determine success based on whether the residual is less than the tolerance
end
