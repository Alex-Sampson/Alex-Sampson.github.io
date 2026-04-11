% task1_3_driver.m  (FIXED: Figure 5 bracket-index plot will always show)
%
% Task 1: f(x) = (x-rho)^d
% Uses your routines (exact signatures):
%   [x, xhist, fhist, dxhist, k] = Mnewton(f, fp, x0, m, tol, maxit)
%   [x, xhist, fhist, dxhist, k] = steffensen(f, x0, tol, maxit)
%   [x, xhist, fhist, dxhist, k] = secant(f, x0, x1, tol, maxit)
%   [x, xhist, fhist, dxhist, k, a, b] = regula_falsi(f, a0, b0, tol, maxit)
%
% Output:
%   3 experiments x (table + figures)  -> 6 figures total
%
% Figure 1: Newton-family convergence histories (|f| vs k) for multiple d
% Figure 2: Iterations vs d (Newton family + Steffensen)
% Figure 3: Secant sensitivity (iterations vs pair index) for multiple d
% Figure 4: Secant residual histories for representative pairs
% Figure 5: Regula Falsi sensitivity (iterations vs bracket index) for odd d
% Figure 6: Regula Falsi residual histories for representative brackets (odd d)

clear; clc; close all;  % Clear the workspace, command window, and close all figures

% ---------------- settings ----------------
rho   = 1.0;  % Set the value of rho, which is a parameter for our functions
tol   = 1e-12;  % Define the tolerance level for convergence
maxit = 60;  % Set the maximum number of iterations allowed

d_list = [1 2 3 4 5 6];  % List of degrees to test

% One-point initial guesses (Newton/Steffensen tables)
x0_list = [rho-0.5, rho-0.1, rho+0.1, rho+0.5];  % Initial guesses for the root

% One x0 used for "rate" plots 
x0_show = rho - 0.5;  % A specific initial guess for plotting purposes

% ------------------------------------------------------------
% IMPORTANT: choose ASYMMETRIC secant pairs so secant doesn't hit rho in 1 step
% (symmetric pairs often make secant land exactly at rho for odd d)
% ------------------------------------------------------------
sec_pairs = [ rho-0.5, rho+0.2;  % Define pairs for the secant method
              rho-1.0, rho+0.5;
              rho-0.2, rho+1.0;
              rho-0.8, rho+0.1;
              rho+0.2, rho+1.0;    % same-side pair
              rho-1.5, rho-0.2 ];  % same-side pair

% Representative pairs to show residual histories
sec_pairs_to_plot = [1 2 3];  % Select pairs to plot their residuals

% Multiplicities to use in secant sensitivity plot
d_sec_list = [2 3 5 6];  % Degrees to test for secant sensitivity

% ------------------------------------------------------------
% choose ASYMMETRIC RF brackets so false-position doesn't hit rho in 1 step
% (symmetric [rho-w, rho+w] gives x=rho immediately for odd d)
% ------------------------------------------------------------
rf_brackets = [ rho-0.5, rho+1.0;  % Define brackets for the Regula Falsi method
                rho-1.0, rho+0.4;
                rho-0.2, rho+0.8;
                rho-2.0, rho+0.5 ];


d_rf_list = [3 5];  % List of odd degrees for Regula Falsi sensitivity


tol_rf   = tol;  % Set the same tolerance for Regula Falsi
maxit_rf = 500;  % Increase max iterations for Regula Falsi to help ensure convergence


tiny_plot = 1e-300;  % A tiny value to avoid log(0) issues in plots

fprintf('\nTASK 1 DRIVER (fixed)\n');  % Print task header
fprintf('rho=%.6g, tol=%.1e, maxit=%d (RF uses maxit_rf=%d)\n\n', rho, tol, maxit, maxit_rf);  % Print settings

% ============================================================
% EXPERIMENT 1: Newton family + Steffensen
%   Table for all d, and Figures 1-2
% ============================================================

fprintf('============================================================\n');  % Print separator
fprintf('EXPERIMENT 1: Newton family + Steffensen\n');  % Print experiment title
fprintf('============================================================\n\n');  % Print separator


for di = 1:length(d_list)  % Loop over each degree in d_list
    d  = d_list(di);  % Get the current degree
    f  = @(x) (x - rho).^d;  % Define the function for the current degree
    fp = @(x) d*(x - rho).^(d-1);  % Define the derivative of the function

    fprintf('--- d=%d ---\n', d);  % Print the current degree
    fprintf('%-18s %-10s %-6s %-6s %-20s %-14s\n', ...  % Print table header
        'Method', 'x0', 'Iter', 'Conv', 'x_end', '|f(x_end)|');
    fprintf('%s\n', repmat('-',1,80));  % Print a separator line

    for j = 1:length(x0_list)  % Loop over each initial guess
        x0 = x0_list(j);  % Get the current initial guess

        % Newton m=1
        [x, ~, fh, ~, k] = Mnewton(f, fp, x0, 1.0, tol, maxit);  % Run Newton's method with m=1
        conv = (~isempty(fh) && abs(fh(end)) < tol);  % Check if it converged
        fprintf('%-18s %-10.4g %-6d %-6d %-20.16f %-14.3e\n', ...  % Print results
            'Newton (m=1)', x0, k, conv, x, abs(f(x)));

        % Newton m=d
        [x, ~, fh, ~, k] = Mnewton(f, fp, x0, double(d), tol, maxit);  % Run Newton's method with m=d
        conv = (~isempty(fh) && abs(fh(end)) < tol);  % Check for convergence
        fprintf('%-18s %-10.4g %-6d %-6d %-20.16f %-14.3e\n', ...  % Print results 
            'Newton (m=d)', x0, k, conv, x, abs(f(x)));

      
        if d >= 2  % Check if d is at least 2
            [x, ~, fh, ~, k] = Mnewton(f, fp, x0, double(d-1), tol, maxit);  % Run Newton's method with m=d-1
            conv = (~isempty(fh) && abs(fh(end)) < tol);  % Check for convergence
            fprintf('%-18s %-10.4g %-6d %-6d %-20.16f %-14.3e\n', ...  % Print results
                'Newton (m=d-1)', x0, k, conv, x, abs(f(x)));

            [x, ~, fh, ~, k] = Mnewton(f, fp, x0, double(d+1), tol, maxit);  % Run Newton's method with m=d+1
            conv = (~isempty(fh) && abs(fh(end)) < tol);  % Check for convergence
            fprintf('%-18s %-10.4g %-6d %-6d %-20.16f %-14.3e\n', ...  % Print results
                'Newton (m=d+1)', x0, k, conv, x, abs(f(x)));
        end

        % Steffensen
        [x, ~, fh, ~, k] = steffensen(f, x0, tol, maxit);  % Run Steffensen's method
        conv = (~isempty(fh) && abs(fh(end)) < tol);  % Check for convergence
        fprintf('%-18s %-10.4g %-6d %-6d %-20.16f %-14.3e\n', ...  % Print results
            'Steffensen', x0, k, conv, x, abs(f(x)));

        fprintf('%s\n', repmat('.',1,80));  % Print a separator line
    end
    fprintf('\n');  % Print a new line for better readability
end


% ---- Figure 1: convergence histories for multiple d (Newton family + Steffensen) ----
d_hist_list = [2 3 5 6];  % Define a list of degrees to analyze convergence

figure(1); clf;  % Create a new figure and clear any existing plots
for p = 1:length(d_hist_list)  % Loop through each degree in the list
    d  = d_hist_list(p);  % Get the current degree
    f  = @(x) (x - rho).^d;  % Define the function for the current degree
    fp = @(x) d*(x - rho).^(d-1);  % Define the derivative of the function

    [~, ~, f1,  ~, ~] = Mnewton(f, fp, x0_show, 1.0,       tol, maxit);  % Run Newton's method with m=1
    [~, ~, fd,  ~, ~] = Mnewton(f, fp, x0_show, double(d), tol, maxit);  % Run Newton's method with m=d

    if d >= 2  % Check if the degree is at least 2
        [~, ~, fm,  ~, ~] = Mnewton(f, fp, x0_show, double(d-1), tol, maxit);  % Run Newton's method with m=d-1
        [~, ~, fp1, ~, ~] = Mnewton(f, fp, x0_show, double(d+1), tol, maxit);  % Run Newton's method with m=d+1
    else
        fm = []; fp1 = [];  % If d is less than 2, set fm and fp1 to empty
    end

    [~, ~, fs, ~, ~] = steffensen(f, x0_show, tol, maxit);  % Run Steffensen's method

    subplot(2,2,p); grid on; hold on;  % Create a subplot for the current degree, enable grid and hold on for multiple plots
    y = max(abs(f1), tiny_plot);   semilogy(0:length(y)-1, y, 'LineWidth', 1.5);  % Plot the absolute values of f1 on a logarithmic scale
    y = max(abs(fd), tiny_plot);   semilogy(0:length(y)-1, y, 'LineWidth', 1.5);  % Plot the absolute values of fd on a logarithmic scale
    if ~isempty(fm),  y = max(abs(fm),  tiny_plot); semilogy(0:length(y)-1, y, 'LineWidth', 1.5); end  % Plot fm if it's not empty
    if ~isempty(fp1), y = max(abs(fp1), tiny_plot); semilogy(0:length(y)-1, y, 'LineWidth', 1.5); end  % Plot fp1 if it's not empty
    y = max(abs(fs), tiny_plot);   semilogy(0:length(y)-1, y, 'LineWidth', 1.5);  % Plot the absolute values of fs on a logarithmic scale

    title(sprintf('d=%d (x0=%.2g)', d, x0_show));  % Set the title for the subplot with the current degree and initial guess
    xlabel('k'); ylabel('|f(x_k)|');  % Label the x-axis and y-axis

    leg = {'Newton m=1','Newton m=d'};  % Initialize the legend with the first two methods
    if ~isempty(fm),  leg{end+1} = 'Newton m=d-1'; end  % Add m=d-1 to the legend if it's not empty
    if ~isempty(fp1), leg{end+1} = 'Newton m=d+1'; end  % Add m=d+1 to the legend if it's not empty
    leg{end+1} = 'Steffensen';  % Add Steffensen's method to the legend
    legend(leg, 'Location','best');  % Display the legend in the best location
    hold off;  % Release the hold on the current subplot
end
sgtitle('Figure 1: Newton-family residual histories (multiple d)');  % Set a super title for the entire figure

% ---- Figure 2: iterations vs d (Newton family + Steffensen) ----
its_m1  = NaN(size(d_list));  % Initialize an array to hold iteration counts for m=1, filled with NaNs
its_md  = NaN(size(d_list));  % Initialize an array for m=d
its_dm1 = NaN(size(d_list));  % Initialize an array for m=d-1
its_st  = NaN(size(d_list));  % Initialize an array for Steffensen's method

x0 = x0_show;  % Set the initial guess for the iterations
for di = 1:length(d_list)  % Loop through each degree in d_list
    d  = d_list(di);  % Get the current degree
    f  = @(x) (x - rho).^d;  % Define the function for the current degree
    fp = @(x) d*(x - rho).^(d-1);  % Define the derivative of the function

    [~, ~, fh, ~, k] = Mnewton(f, fp, x0, 1.0, tol, maxit);  % Run Newton's method with m=1
    if ~isempty(fh) && abs(fh(end)) < tol, its_m1(di) = k; end  % Store the iteration count if it converged

    [~, ~, fh, ~, k] = Mnewton(f, fp, x0, double(d), tol, maxit);  % Run Newton's method with m=d
    if ~isempty(fh) && abs(fh(end)) < tol, its_md(di) = k; end  % Store the iteration count if it converged

    if d >= 2  % Check if the degree is at least 2
        [~, ~, fh, ~, k] = Mnewton(f, fp, x0, double(d-1), tol, maxit);  % Run Newton's method with m=d-1
        if ~isempty(fh) && abs(fh(end)) < tol, its_dm1(di) = k; end  % Store the iteration count if it converged
    end

    [~, ~, fh, ~, k] = steffensen(f, x0, tol, maxit);  % Run Steffensen's method
    if ~isempty(fh) && abs(fh(end)) < tol, its_st(di) = k; end  % Store the iteration count if it converged
end

figure(2); clf; grid on; hold on;  % Create a new figure for the iterations vs d, clear it, enable grid and hold on for multiple plots
plot(d_list, its_m1,  '-o', 'LineWidth', 1.5);  % Plot the iterations for m=1 against the degrees
plot(d_list, its_md,  '-o', 'LineWidth', 1.5);  % Plot the iterations for m=d against the degrees
plot(d_list, its_dm1, '-o', 'LineWidth', 1.5);  % Plot the iterations for m=d-1 against the degrees
plot(d_list, its_st,  '-o', 'LineWidth', 1.5);  % Plot the iterations for Steffensen's method against the degrees
xticks(d_list);  % Set the x-ticks to the degrees in d_list
xlabel('d'); ylabel('iterations (NaN = failed)');  % Label the x-axis and y-axis
title(sprintf('Figure 2: Iterations vs d (x0=%.2g)', x0));  % Set the title for the figure with the initial guess
legend('Newton m=1','Newton m=d','Newton m=d-1','Steffensen','Location','best');  % Create a legend for the different methods
ylim([0, maxit]);  % Set the y-axis limits to be between 0 and the maximum iterations
hold off;  % Release the hold on the current figure

% ============================================================
% EXPERIMENT 2: Secant sensitivity
%   Figures 3-4
% ============================================================

fprintf('============================================================\n');  % Print a separator line for clarity
fprintf('EXPERIMENT 2: Secant sensitivity\n');  % Indicate the start of Experiment 2
fprintf('============================================================\n\n');  % Print another separator line

% ---- Figure 3: iterations vs pair index for multiple d ----
figure(3); clf; grid on; hold on;  % Create a new figure, clear it, and set up the grid for plotting
for di = 1:length(d_sec_list)  % Loop through each degree in the list of secant degrees
    d = d_sec_list(di);  % Get the current degree
    f = @(x) (x - rho).^d;  % Define the function for the current degree

    its = NaN(size(sec_pairs,1),1);  % Initialize an array to hold iteration counts, filled with NaNs
    for j = 1:size(sec_pairs,1)  % Loop through each pair of starting points
        [~, ~, fh, ~, k] = secant(f, sec_pairs(j,1), sec_pairs(j,2), tol, maxit);  % Run the secant method on the current pair
        if ~isempty(fh) && abs(fh(end)) < tol  % Check if the history is not empty and if the last value is within tolerance
            its(j) = k;  % Store the number of iterations if it converged
        else
            its(j) = NaN;  % Otherwise, mark it as NaN to indicate failure
        end
    end
    plot(1:size(sec_pairs,1), its, '-o', 'LineWidth', 1.5);  % Plot the iterations against the pair index
end
xlabel('pair index'); ylabel('iterations (NaN = failed)');  % Label the axes
title('Figure 3: Secant sensitivity (iterations vs starting pair)');  % Set the title for the figure
legend(arrayfun(@(dd) sprintf('d=%d',dd), d_sec_list, 'UniformOutput', false), ...  % Create a legend for the different degrees
       'Location','best');  % Place the legend in the best location
ylim([0, maxit]);  % Set the y-axis limits to be between 0 and the maximum iterations

% --- reserve space on the right for the pair index key ---
set(gca,'Position',[0.08 0.11 0.58 0.8])  % Adjust the position of the axes to make room for annotations

% --- build text describing which pair corresponds to each index ---
pair_text = cell(size(sec_pairs,1),1);  % Create a cell array to hold the text for each pair
for j = 1:size(sec_pairs,1)  % Loop through each pair
    pair_text{j} = sprintf('%d: (%.2g, %.2g)', j, sec_pairs(j,1), sec_pairs(j,2));  % Format the text to show the index and the pair values
end

% --- place annotation outside the plotting area ---
annotation('textbox', [0.70 0.15 0.28 0.7], ...  % Create a textbox annotation
    'String', pair_text, ...  % Set the content of the textbox to the pair text
    'FitBoxToText', 'on', ...  % Adjust the box to fit the text
    'BackgroundColor', 'white', ...  % Set the background color of the textbox to white
    'EdgeColor', 'black', ...  % Set the edge color of the textbox to black
    'FontSize', 9);  % Set the font size of the text in the textbox

hold off;  % Release the hold on the current figure

% ---- Figure 4: residual histories for representative pairs (two d values) ----
dA = 2; dB = 5;  % Define two degrees for comparison
figure(4); clf;  % Create a new figure for the residual histories and clear it

% Top: d=2
subplot(2,1,1); grid on; hold on;  % Create the first subplot and set up the grid
f = @(x) (x - rho).^dA;  % Define the function for degree dA
for t = 1:length(sec_pairs_to_plot)  % Loop through the pairs to plot
    j = sec_pairs_to_plot(t);  % Get the current pair index
    [~, ~, fh, ~, ~] = secant(f, sec_pairs(j,1), sec_pairs(j,2), tol, maxit);  % Run the secant method on the current pair
    if ~isempty(fh)  % Check if the history is not empty
        y = max(abs(fh), tiny_plot);  % Get the maximum absolute value of the history, ensuring it's above a tiny threshold
        semilogy(0:length(y)-1, y, 'LineWidth', 1.5);  % Plot the residuals on a logarithmic scale
    end
end
title(sprintf('Figure 4A: Secant residual histories (d=%d)', dA));  % Set the title for the first subplot
ylabel('|f(x_k)|');  % Label the y-axis
legend(arrayfun(@(j) sprintf('pair %d',j), sec_pairs_to_plot, 'UniformOutput', false), ...  % Create a legend for the pairs being plotted
       'Location','best');  % Place the legend in the best location
hold off;  % Release the hold on the current subplot

% Bottom: d=5
subplot(2,1,2); grid on; hold on;  % Create the second subplot and set up the grid
f = @(x) (x - rho).^dB;  % Define the function for degree dB
for t = 1:length(sec_pairs_to_plot)  % Loop through the pairs to plot
    j = sec_pairs_to_plot(t);  % Get the current pair index
    [~, ~, fh, ~, ~] = secant(f, sec_pairs(j,1), sec_pairs(j,2), tol, maxit);  % Run the secant method on the current pair
    if ~isempty(fh)  % Check if the history is not empty
        y = max(abs(fh), tiny_plot);  % Get the maximum absolute value of the history, ensuring it's above a tiny threshold
        semilogy(0:length(y)-1, y, 'LineWidth', 1.5);  % Plot the residuals on a logarithmic scale
    end
end
title(sprintf('Figure 4B: Secant residual histories (d=%d)', dB));  % Set the title for the second subplot
xlabel('k'); ylabel('|f(x_k)|');  % Label the x and y axes
legend(arrayfun(@(j) sprintf('pair %d',j), sec_pairs_to_plot, 'UniformOutput', false), ...  % Create a legend for the pairs being plotted
       'Location','best');  % Place the legend in the best location
hold off;  % Release the hold on the current subplot

% ============================================================
% EXPERIMENT 3: Regula Falsi sensitivity (odd d only)
%   Figures 5-6
% ============================================================

fprintf('============================================================\n');  % Print a separator line for clarity
fprintf('EXPERIMENT 3: Regula Falsi sensitivity (odd d only)\n');  % Indicate the start of Experiment 3
fprintf('============================================================\n\n');  % Print another separator line

% ---- Figure 5: Regula Falsi final residual vs bracket index (odd d only) ----
figure(5); clf; grid on; hold on;  % Create a new figure, clear it, and set up the grid for plotting

for di = 1:length(d_rf_list)  % Loop through each degree in the list of odd degrees for Regula Falsi
    d = d_rf_list(di);  % Get the current degree
    f = @(x) (x - rho).^d;  % Define the function for the current degree

    final_res = NaN(size(rf_brackets,1),1);  % Initialize an array to hold final residuals, filled with NaNs

    for j = 1:size(rf_brackets,1)  % Loop through each bracket in the list
        a0 = rf_brackets(j,1);  % Get the lower bound of the current bracket
        b0 = rf_brackets(j,2);  % Get the upper bound of the current bracket

        if f(a0)*f(b0) <= 0  % Check if the function values at the bracket ends have opposite signs
            [~, ~, fh, ~, ~] = regula_falsi(f, a0, b0, tol_rf, maxit_rf);  % Run the Regula Falsi method on the current bracket
            if ~isempty(fh)  % Check if the history of function values is not empty
                [x_end, xh, fh, ~, ~] = regula_falsi(f, a0, b0, tol_rf, maxit_rf);  % Run the method again to get the final results
                    if ~isempty(xh)  % Check if the history of x values is not empty
                        final_res(j) = abs(xh(end) - rho);   % Calculate the distance to the root and store it
                    end

            end
        end
    end

    semilogy(1:size(rf_brackets,1), final_res, '-o', 'LineWidth', 1.5);  % Plot the final residuals on a logarithmic scale
end

xlabel('bracket index');  % Label the x-axis
ylabel(sprintf('|x_{final}-\\rho| after %d iterations', maxit_rf));  % Label the y-axis with the final residual description
title('Figure 5: Regula Falsi sensitivity to bracket choice (odd d only)');  % Set the title for the figure
legend(arrayfun(@(dd) sprintf('d=%d',dd), d_rf_list, 'UniformOutput', false), ...  % Create a legend for the different degrees
       'Location','best');  % Place the legend in the best location
hold off;  % Release the hold on the current figure

% ---- Figure 6: RF residual histories for two brackets (odd d, choose d=5) ----
d = 5;  % Set the degree for the residual history plots
f = @(x) (x - rho).^d;  % Define the function for degree d

br1 = 1;  % First bracket in the list to plot
br2 = 4;  % Last bracket in the list to plot

a1 = rf_brackets(br1,1); b1 = rf_brackets(br1,2);  % Get the bounds for the first bracket
a2 = rf_brackets(br2,1); b2 = rf_brackets(br2,2);  % Get the bounds for the second bracket

figure(6); clf;  % Create a new figure for the residual histories and clear it

subplot(2,1,1); grid on; hold on;  % Create the first subplot and set up the grid
if f(a1)*f(b1) <= 0  % Check if the function values at the first bracket ends have opposite signs
    [~, ~, fh, ~, ~, ~, ~] = regula_falsi(f, a1, b1, tol_rf, maxit_rf);  % Run the Regula Falsi method on the first bracket
    if ~isempty(fh)  % Check if the history of function values is not empty
        y = max(abs(fh), tiny_plot);  % Get the maximum absolute value of the history, ensuring it's above a tiny threshold
        semilogy(0:length(y)-1, y, 'LineWidth', 1.5);  % Plot the residuals on a logarithmic scale
    end
end
title(sprintf('Figure 6A: RF residual history (d=%d, bracket %d)', d, br1));  % Set the title for the first subplot
ylabel('|f(x_k)|');  % Label the y-axis
hold off;  % Release the hold on the current subplot

subplot(2,1,2); grid on; hold on;  % Create the second subplot and set up the grid
if f(a2)*f(b2) <= 0  % Check if the function values at the second bracket ends have opposite signs
    [~, ~, fh, ~, ~, ~, ~] = regula_falsi(f, a2, b2, tol_rf, maxit_rf);  % Run the Regula Falsi method on the second bracket
    if ~isempty(fh)  % Check if the history of function values is not empty
        y = max(abs(fh), tiny_plot);  % Get the maximum absolute value of the history, ensuring it's above a tiny threshold
        semilogy(0:length(y)-1, y, 'LineWidth', 1.5);  % Plot the residuals on a logarithmic scale
    end
end
title(sprintf('Figure 6B: RF residual history (d=%d, bracket %d)', d, br2));  % Set the title for the second subplot
xlabel('k'); ylabel('|f(x_k)|');  % Label the x and y axes
hold off;  % Release the hold on the current subplot

fprintf('DONE. Produced Figures 1–6.\n');  % Print a completion message
