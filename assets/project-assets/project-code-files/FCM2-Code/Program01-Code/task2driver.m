% task2_driver.m
%
% Task 2 driver for Programming Assignment 1 (FCM2 Spring 2026)
%
% UPDATE (Fix): Some plots were labeling "(failures shown at 61)" even when
% there were no failures. This version sets that label ONLY if failures occur.
%
% Fixes included:
%   (A) Excursion size = max_{k>=0} |x_k| INCLUDING x0.
%   (B) Task 2.2 sigma axis uses sigma = logspace(-16,16,...) so log10(sigma) in [-16,16].
%   (C) Task 2.3 distance-to-nearest-root histories: full histories + proper log plotting,
%       and fixed-step Newton so tiny rho1 doesn't instantly satisfy |f(x0)|<tol at k=0.
%   (D) NEW: only show "(failures shown at maxit+1)" when failures actually exist.
%
% Requires in same folder:
%   newton.m, steffensen.m, secant.m, regula_falsi.m
%
% ------------------------------------------------------------

clear; clc; % Clear the workspace and command window for a fresh start

%% user parameters
rho   = 1.0;     % Set the value of rho, which is used for defining roots in the function
tol   = 1e-12;   % Define the tolerance level for convergence
maxit = 60;      % Set the maximum number of iterations allowed for the methods

% plotting range for x0 sweeps (Task 2.1(3))
x0_min = -2*rho; % Minimum value for initial guess x0, set to -2 times rho
x0_max =  2*rho; % Maximum value for initial guess x0, set to +2 times rho
Nx0    = 3001;    % Number of points in the x0 grid for sweeping

%% Task 2.1 setup: f(x)=x(x-rho)(x+rho)=x^3-rho^2 x
f  = @(x) x.^3 - (rho^2)*x; % Define the function f(x) using an anonymous function
fp = @(x) 3*x.^2 - rho^2;    % Define the derivative of f(x) using an anonymous function

alpha = rho/sqrt(3);   % Calculate alpha, where the derivative is zero at +/-alpha
xi    = rho/sqrt(5);   % Calculate xi, the symmetric 2-cycle point for this function

roots_21 = [-rho, 0, rho]; % Define the roots of the function in an array

fprintf('============================================================\n'); % Print a separator line
fprintf('TASK 2.1: f(x)=x(x-rho)(x+rho), rho=%.4g\n', rho); % Display the task and value of rho
fprintf('roots: -rho, 0, +rho  =  %.4g, %.4g, %.4g\n', -rho, 0, rho); % Show the roots
fprintf('alpha = rho/sqrt(3) = %.16g (f''=0 at +/-alpha)\n', alpha); % Print the value of alpha
fprintf('xi    = rho/sqrt(5) = %.16g (symmetric 2-cycle point)\n', xi); % Print the value of xi
fprintf('tol=%.1e, maxit=%d\n', tol, maxit); % Show the tolerance and maximum iterations
fprintf('============================================================\n\n'); % Print another separator line



%% ============================================================
% local helper functions
% ============================================================

function [ok, root_hit] = classify_convergence(xend, f, tol, roots_list)
    ok = is_converged(xend, f, tol); % Check if the final value xend has converged based on the function f and tolerance tol
    root_hit = NaN; % Initialize root_hit as NaN, will store the root if convergence is successful
    if ok % If the convergence check passed
        [~, idx] = min(abs(xend - roots_list)); % Find the index of the closest root to xend
        root_hit = roots_list(idx); % Set root_hit to the closest root found
    end
end

function ok = is_converged(xend, f, tol)
    if ~isfinite(xend) % Check if xend is a finite number
        ok = false; return; % If not, set ok to false and exit
    end
    fx = f(xend); % Evaluate the function f at xend
    ok = isfinite(fx) && abs(fx) < tol; % Check if fx is finite and within the tolerance
end

function k = run_count_only(newton_fun, f, fp, x0, tol, maxit)
    [xend, xhist, fhist, dxhist, k] = newton_fun(f, fp, x0, tol, maxit); % Run the Newton method and get results
    if ~isfinite(xend) || ~isfinite(f(xend)) || abs(f(xend)) >= tol % Check if the results are valid and within tolerance
        k = NaN; % If not, set k to NaN to indicate failure
    end
end

function k = run_count_only_secant(f, x0, x1, tol, maxit)
    [xend, xhist, fhist, dxhist, k] = secant(f, x0, x1, tol, maxit); % Run the secant method and get results
    if ~isfinite(xend) || ~isfinite(f(xend)) || abs(f(xend)) >= tol % Check if the results are valid and within tolerance
        k = NaN; % If not, set k to NaN to indicate failure
    end
end

function y = replace_fail(y, maxit)
    y(~isfinite(y)) = maxit + 1; % Replace non-finite values in y with maxit + 1 to indicate failure
end

function xhist = newton_fixed_steps(f, fp, x0, K)
    % fixed-step Newton (no convergence stopping), with basic blow-up checks
    tiny = 1e-14; % Define a small threshold for checking the derivative
    x = x0; % Initialize x with the starting point x0
    xhist = x; % Start the history of x values with the initial value

    for k = 1:K % Loop for a fixed number of steps K
        fx  = f(x); % Evaluate the function f at the current x
        dfx = fp(x); % Evaluate the derivative of f at the current x

        if ~isfinite(x) || ~isfinite(fx) || ~isfinite(dfx) % Check if x, fx, or dfx are finite
            break; % If any are not finite, exit the loop
        end
        if abs(dfx) < tiny % Check if the derivative is too small
            break; % If it is, exit the loop to avoid division by zero
        end

        dx = -fx/dfx; % Calculate the Newton step
        x  = x + dx; % Update x with the new value
        xhist = [xhist; x]; % Append the new x value to the history
    end
end



%% ------------------------------------------------------------
% Task 2.1(1): demo runs
% ------------------------------------------------------------
x0_demo = [ 2*rho,  0.8*rho,  0.2*rho,  0.05*rho, ... % Define a set of initial guesses for demo runs
           -0.05*rho,-0.2*rho,-0.8*rho,-2*rho];

fprintf('Task 2.1(1) demo runs:\n'); % Print header for demo runs
fprintf('   x0          k     x_end            f(x_end)          root_hit\n'); % Print column headers
fprintf('---------------------------------------------------------------\n'); % Print a separator line
for j = 1:numel(x0_demo) % Loop through each initial guess in the demo array
    x0 = x0_demo(j); % Set the current initial guess
    [xend, xhist, fhist, dxhist, k] = newton(f, fp, x0, tol, maxit); % Call the Newton method
    [ok, root_hit] = classify_convergence(xend, f, tol, roots_21); % Check if the method converged and which root was hit
    if ok % If the method converged successfully
        fprintf('%+11.4e  %3d  %+ .16e  %+ .2e   %+g\n', x0, k, xend, f(xend), root_hit); % Print results
    else % If the method failed to converge
        fprintf('%+11.4e  %3d  %+ .16e  %+ .2e   %s\n', x0, k, xend, f(xend), 'FAIL'); % Indicate failure
    end
end
fprintf('\n'); % Print a newline for better readability

%% ------------------------------------------------------------
% Task 2.1(2): cycling behavior near xi
% ------------------------------------------------------------
deltas = [0, 1e-16, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6]; % Define small perturbations around xi

figure(1); clf; hold on; grid on; % Create a new figure and set up for plotting
title('Task 2.1(2): Newton iterates near the cycling point xi'); % Set the title of the plot
xlabel('iteration k'); % Label the x-axis
ylabel('x_k'); % Label the y-axis

for j = 1:numel(deltas) % Loop through each delta value
    x0 = xi + deltas(j); % Set the initial guess by adding delta to xi
    [xend, xhist, fhist, dxhist, k] = newton(f, fp, x0, tol, maxit); % Call the Newton method
    plot(0:(numel(xhist)-1), xhist, '-o', 'DisplayName', sprintf('x0 = xi %+g', deltas(j))); % Plot the history of x values
end

yline(+xi,'--','DisplayName','+xi'); % Draw a horizontal line at +xi
yline(-xi,'--','DisplayName','-xi'); % Draw a horizontal line at -xi
legend('Location','northeast'); % Add a legend to the plot
hold off; % Release the hold on the current figure

%% ------------------------------------------------------------
% Task 2.1(3): sweep x0 and plot:
%   (i) which root is hit
%   (ii) iteration count
%   (iii) excursion size E(x0)=max_k |x_k| including k=0
% ------------------------------------------------------------
x0_grid = linspace(x0_min, x0_max, Nx0); % Create a grid of initial guesses between x0_min and x0_max

root_hit = NaN(size(x0_grid)); % Initialize an array to store which root was hit
iters    = NaN(size(x0_grid)); % Initialize an array to store the number of iterations
excurs   = NaN(size(x0_grid)); % Initialize an array to store the excursion sizes

n_fail = 0; % Initialize a counter for failed runs

for i = 1:numel(x0_grid) % Loop through each initial guess in the grid
    x0 = x0_grid(i); % Set the current initial guess
    [xend, xhist, fhist, dxhist, k] = newton(f, fp, x0, tol, maxit); % Call the Newton method

    [ok, r] = classify_convergence(xend, f, tol, roots_21); % Check if the method converged and which root was hit

    if ok % If the method converged successfully
        root_hit(i) = r; % Store the root that was hit
        iters(i)    = k; % Store the number of iterations taken
        % FIX (A): include x0 by taking max over FULL xhist
        excurs(i)   = max(abs(xhist)); % Calculate the excursion size including the initial guess
    else % If the method failed to converge
        n_fail = n_fail + 1; % Increment the failure counter
    end
end

figure(2); clf; hold on; grid on; % Create a new figure for plotting and set up the grid
title('Task 2.1(3): Newton root hit vs x0'); % Set the title of the plot
xlabel('initial guess x0'); % Label the x-axis
ylabel('root hit (approx: -rho, 0, +rho)'); % Label the y-axis
plot(x0_grid, root_hit, '.', 'MarkerSize', 8); % Plot the roots hit for each initial guess
xline(-rho,'--'); xline(0,'--'); xline(+rho,'--'); % Draw vertical lines at the roots
xline(-alpha,':'); xline(+alpha,':'); % Draw dashed lines at alpha
xline(-xi,'-.');   xline(+xi,'-.'); % Draw dash-dot lines at xi
yticks([-rho 0 rho]); % Set the y-ticks to the roots
hold off; % Release the hold on the current figure

figure(3); clf; hold on; grid on; % Create a new figure for plotting and set up the grid
title('Task 2.1(3): Newton iteration count vs x0'); % Set the title of the plot
xlabel('initial guess x0'); % Label the x-axis
ylabel('Newton iterations (converged runs)'); % Label the y-axis
plot(x0_grid, iters, '.', 'MarkerSize', 8); % Plot the iteration counts for each initial guess
xline(-rho,'--'); xline(0,'--'); xline(+rho,'--'); % Draw vertical lines at the roots
xline(-alpha,':'); xline(+alpha,':'); % Draw dashed lines at alpha
xline(-xi,'-.');   xline(+xi,'-.'); % Draw dash-dot lines at xi
hold off; % Release the hold on the current figure

figure(4); clf; hold on; grid on; % Create a new figure for plotting and set up the grid
title('Task 2.1(3): Newton excursion size vs x0'); % Set the title of the plot
xlabel('initial guess x0'); % Label the x-axis
ylabel('max_k |x_k| over iterations (converged runs)'); % Label the y-axis
set(gca,'YScale','log'); % Set the y-axis to a logarithmic scale
plot(x0_grid, excurs, '.', 'MarkerSize', 8); % Plot the excursion sizes for each initial guess
xline(-rho,'--'); xline(0,'--'); xline(+rho,'--'); % Draw vertical lines at the roots
xline(-alpha,':'); xline(+alpha,':'); % Draw dashed lines at alpha
xline(-xi,'-.');   xline(+xi,'-.'); % Draw dash-dot lines at xi
hold off; % Release the hold on the current figure

fprintf('Task 2.1(3): %d/%d x0 values failed (nonconvergent/early abort).\n\n', n_fail, numel(x0_grid)); % Print the number of failed runs

%% ------------------------------------------------------------
% Task 2.2: scaling: f~(x)=sigma*f(x)
% Compare Newton/Steffensen/Secant/Regula Falsi as sigma varies
% ------------------------------------------------------------
fprintf('============================================================\n'); % Print a separator line for clarity
fprintf('TASK 2.2: Scaling  f~(x) = sigma*f(x)\n'); % Announce the task being performed
fprintf('============================================================\n'); % Print another separator line

% FIX (B): sensible sigma range for log10 axis
sigma_list = logspace(-16, 16, 9);   % Create a list of sigma values on a logarithmic scale from 10^-16 to 10^16

% initial conditions
x0_newt = 0.8*rho;          % Set the initial guess for Newton's method close to the root
x0_sec0 = -1.2*rho;         % Set the first initial guess for the secant method, away from the root
x0_sec1 =  0.8*rho;         % Set the second initial guess for the secant method, close to the root
a_rf    = -1.5*rho;         % Define the lower bound for the Regula Falsi method
b_rf    =  0.2*rho;         % Define the upper bound for the Regula Falsi method

iters_newt = NaN(size(sigma_list)); % Preallocate an array to store iteration counts for Newton's method
iters_stef = NaN(size(sigma_list)); % Preallocate an array for Steffensen's method
iters_sec  = NaN(size(sigma_list)); % Preallocate an array for the Secant method
iters_rf   = NaN(size(sigma_list)); % Preallocate an array for the Regula Falsi method

for j = 1:numel(sigma_list) % Loop through each sigma value
    sigma = sigma_list(j); % Get the current sigma value from the list
    fs  = @(x) sigma * f(x); % Define the scaled function for Newton's method
    fps = @(x) sigma * fp(x); % Define the scaled derivative for Newton's method

    % Newton
    [xend, xhist, fhist, dxhist, k] = newton(fs, fps, x0_newt, tol, maxit); % Call Newton's method with the scaled function
    iters_newt(j) = k; % Store the number of iterations taken
    if ~is_converged(xend, fs, tol), iters_newt(j) = NaN; end % Check if the method converged, if not set to NaN

    % Steffensen
    [xend, xhist, fhist, dxhist, k] = steffensen(fs, x0_newt, tol, maxit); % Call Steffensen's method with the scaled function
    iters_stef(j) = k; % Store the number of iterations taken
    if ~is_converged(xend, fs, tol), iters_stef(j) = NaN; end % Check for convergence, set to NaN if not

    % Secant
    [xend, xhist, fhist, dxhist, k] = secant(fs, x0_sec0, x0_sec1, tol, maxit); % Call the Secant method with the scaled function
    iters_sec(j) = k; % Store the number of iterations taken
    if ~is_converged(xend, fs, tol), iters_sec(j) = NaN; end % Check for convergence, set to NaN if not

    % Regula Falsi
    [xend, xhist, fhist, dxhist, k, aa, bb] = regula_falsi(fs, a_rf, b_rf, tol, maxit); % Call the Regula Falsi method with the scaled function
    iters_rf(j) = k; % Store the number of iterations taken
    if ~is_converged(xend, fs, tol), iters_rf(j) = NaN; end % Check for convergence, set to NaN if not
end

hasFail_22 = any(~isfinite(iters_newt)) || any(~isfinite(iters_stef)) || ... % Check if any method failed to converge
             any(~isfinite(iters_sec))  || any(~isfinite(iters_rf)); % Check for failures in all methods

xplot = log10(sigma_list); % Prepare the x-axis values for plotting by taking the log10 of sigma

if hasFail_22 % If there were any failures in the methods
    yN  = iters_newt; yN(~isfinite(yN)) = maxit+1; % Replace NaNs with a value indicating failure for Newton
    ySt = iters_stef; ySt(~isfinite(ySt)) = maxit+1; % Replace NaNs for Steffensen
    ySe = iters_sec;  ySe(~isfinite(ySe)) = maxit+1; % Replace NaNs for Secant
    yRf = iters_rf;   yRf(~isfinite(yRf)) = maxit+1; % Replace NaNs for Regula Falsi
else
    yN  = iters_newt; % Use the iteration counts directly for Newton
    ySt = iters_stef; % Use the iteration counts directly for Steffensen
    ySe = iters_sec;  % Use the iteration counts directly for Secant
    yRf = iters_rf;   % Use the iteration counts directly for Regula Falsi
end

figure(5); clf; hold on; grid on; % Create a new figure for plotting and set up the grid
title('Task 2.2: Iterations vs scaling sigma'); % Set the title of the plot
xlabel('log10(sigma)'); % Label the x-axis
if hasFail_22 % If there were any failures
    ylabel(sprintf('iterations (failures shown at %d)', maxit+1)); % Label the y-axis indicating failures
else
    ylabel('iterations'); % Label the y-axis normally
end
plot(xplot, yN,  '-o', 'DisplayName','Newton'); % Plot the iterations for Newton's method
plot(xplot, ySt, '-o', 'DisplayName','Steffensen'); % Plot the iterations for Steffensen's method
plot(xplot, ySe, '-o', 'DisplayName','Secant'); % Plot the iterations for the Secant method
plot(xplot, yRf, '-o', 'DisplayName','Regula Falsi'); % Plot the iterations for the Regula Falsi method
legend('Location','northeast'); % Add a legend to the plot
hold off; % Release the hold on the current figure

%% ------------------------------------------------------------
% Task 2.3: Root coalescing for fhat(x)=x(x-rho1)(x-rho2), 0<rho1<rho2
%  (1) rho1 -> 0 with rho2 fixed
%  (2) rho1,rho2 -> 0 together (rho2 = 2*rho1)
%  plus distance-to-nearest-root history plot (fixed steps)
% ------------------------------------------------------------
fprintf('============================================================\n'); % Print a separator line for clarity
fprintf('TASK 2.3: Root coalescing  fhat(x)=x(x-rho1)(x-rho2)\n'); % Announce the task being performed
fprintf('============================================================\n'); % Print another separator line

% expanded derivative:
% f = x^3 - (r1+r2)x^2 + (r1*r2)x
% f' = 3x^2 - 2(r1+r2)x + (r1*r2)

% ---- (1) rho1 -> 0 with rho2 fixed
rho2_fixed = 1.0; % Set a fixed value for rho2
rho1_list  = logspace(-12, -2, 11); % Create a list of rho1 values spaced logarithmically

iters_near0    = NaN(size(rho1_list)); % Initialize an array to store iteration counts for near 0
iters_between  = NaN(size(rho1_list)); % Initialize an array for iteration counts between roots
iters_near_r2  = NaN(size(rho1_list)); % Initialize an array for iteration counts near rho2
iters_global   = NaN(size(rho1_list)); % Initialize an array for global iteration counts
iters_sec_rep  = NaN(size(rho1_list)); % Initialize an array for Secant method iteration counts

for j = 1:numel(rho1_list) % Loop over each value of rho1
    r1 = rho1_list(j); % Get the current value of rho1
    r2 = rho2_fixed; % Use the fixed value of rho2

    fhat  = @(x) x.*(x-r1).*(x-r2); % Define the function fhat based on current r1 and r2
    fhatp = @(x) 3*x.^2 - 2*(r1+r2)*x + (r1*r2); % Define the derivative of fhat

    x0_near0   = 0.2*r1; % Set an initial guess close to 0
    x0_between = 0.5*(r1+r2); % Set an initial guess between r1 and r2
    x0_near_r2 = r2 - 0.2; % Set an initial guess close to r2
    x0_global  = 2.0*r2; % Set a global initial guess far from the roots

    iters_near0(j)   = run_count_only(@newton, fhat, fhatp, x0_near0,   tol, maxit); % Count iterations for near 0
    iters_between(j) = run_count_only(@newton, fhat, fhatp, x0_between, tol, maxit); % Count iterations for between roots
    iters_near_r2(j) = run_count_only(@newton, fhat, fhatp, x0_near_r2, tol, maxit); % Count iterations for near r2
    iters_global(j)  = run_count_only(@newton, fhat, fhatp, x0_global,  tol, maxit); % Count iterations for global guess

    iters_sec_rep(j) = run_count_only_secant(fhat, x0_between, x0_global, tol, maxit); % Count iterations for Secant method
end

hasFail_231 = any(~isfinite(iters_near0)) || any(~isfinite(iters_between)) || ... % Check if any method failed for near 0 or between
              any(~isfinite(iters_near_r2)) || any(~isfinite(iters_global)) || ... % Check for failures near r2 and global
              any(~isfinite(iters_sec_rep)); % Check for failures in Secant method

if hasFail_231 % If any method failed
    yA = replace_fail(iters_near0,   maxit); % Replace failures for near 0 with a maxit indicator
    yB = replace_fail(iters_between, maxit); % Replace failures for between with a maxit indicator
    yC = replace_fail(iters_near_r2, maxit); % Replace failures for near r2 with a maxit indicator
    yD = replace_fail(iters_global,  maxit); % Replace failures for global with a maxit indicator
    yE = replace_fail(iters_sec_rep, maxit); % Replace failures for Secant with a maxit indicator
else
    yA = iters_near0; yB = iters_between; yC = iters_near_r2; yD = iters_global; yE = iters_sec_rep; % Use original counts if no failures
end

figure(6); clf; hold on; grid on; % Create a new figure for plotting and set up the grid
title('Task 2.3(1): Iterations as rho1 -> 0 with rho2 fixed'); % Set the title of the plot
xlabel('rho1 (log scale), rho2 fixed'); % Label the x-axis
if hasFail_231 % If there were any failures
    ylabel(sprintf('iterations (failures shown at %d)', maxit+1)); % Label the y-axis indicating failures
else
    ylabel('iterations'); % Label the y-axis normally
end
set(gca,'XScale','log'); % Set the x-axis to logarithmic scale
plot(rho1_list, yA, '-o', 'DisplayName','Newton: near0'); % Plot iterations for Newton's method near 0
plot(rho1_list, yB, '-o', 'DisplayName','Newton: between'); % Plot iterations for Newton's method between roots
plot(rho1_list, yC, '-o', 'DisplayName','Newton: near rho2'); % Plot iterations for Newton's method near rho2
plot(rho1_list, yD, '-o', 'DisplayName','Newton: global'); % Plot iterations for Newton's method global
plot(rho1_list, yE, '-o', 'DisplayName','Secant (rep pair)'); % Plot iterations for Secant method
legend('Location','northeast'); % Add a legend to the plot
hold off; % Release the hold on the current figure

% ---- (2) rho1,rho2 -> 0 together, rho2 = 2*rho1
rho1_list2 = logspace(-12, -2, 11); % Create a second list of rho1 values spaced logarithmically

iters_near0_2    = NaN(size(rho1_list2)); % Initialize an array for iteration counts near 0 for the second case
iters_between_2  = NaN(size(rho1_list2)); % Initialize an array for iteration counts between roots for the second case
iters_near_r2_2  = NaN(size(rho1_list2)); % Initialize an array for iteration counts near rho2 for the second case
iters_global_2   = NaN(size(rho1_list2)); % Initialize an array for global iteration counts for the second case
iters_sec_rep_2  = NaN(size(rho1_list2)); % Initialize an array for Secant method iteration counts for the second case

for j = 1:numel(rho1_list2) % Loop over each value of rho1 in the second case
    r1 = rho1_list2(j); % Get the current value of rho1
    r2 = 2*r1; % Set rho2 to be twice the current value of rho1

    fhat  = @(x) x.*(x-r1).*(x-r2); % Define the function fhat based on current r1 and r2
    fhatp = @(x) 3*x.^2 - 2*(r1+r2)*x + (r1*r2); % Define the derivative of fhat

    x0_near0   = 0.2*r1; % Set an initial guess close to 0
    x0_between = 0.5*(r1+r2); % Set an initial guess between r1 and r2
    x0_near_r2 = 0.8*r2; % Set an initial guess close to r2
    x0_global  = 2.0*r2; % Set a global initial guess far from the roots

    iters_near0_2(j)   = run_count_only(@newton, fhat, fhatp, x0_near0,   tol, maxit); % Count iterations for near 0 in the second case
    iters_between_2(j) = run_count_only(@newton, fhat, fhatp, x0_between, tol, maxit); % Count iterations for between roots in the second case
    iters_near_r2_2(j) = run_count_only(@newton, fhat, fhatp, x0_near_r2, tol, maxit); % Count iterations for near r2 in the second case
    iters_global_2(j)  = run_count_only(@newton, fhat, fhatp, x0_global,  tol, maxit); % Count iterations for global guess in the second case

    iters_sec_rep_2(j) = run_count_only_secant(fhat, x0_between, x0_global, tol, maxit); % Count iterations for Secant method in the second case
end

hasFail_232 = any(~isfinite(iters_near0_2)) || any(~isfinite(iters_between_2)) || ... % Check if any method failed for near 0 or between in the second case
              any(~isfinite(iters_near_r2_2)) || any(~isfinite(iters_global_2)) || ... % Check for failures near r2 and global in the second case
              any(~isfinite(iters_sec_rep_2)); % Check for failures in Secant method in the second case

if hasFail_232 % If any method failed in the second case
    yA2 = replace_fail(iters_near0_2,   maxit); % Replace failures for near 0 with a maxit indicator
    yB2 = replace_fail(iters_between_2, maxit); % Replace failures for between with a maxit indicator
    yC2 = replace_fail(iters_near_r2_2, maxit); % Replace failures for near r2 with a maxit indicator
    yD2 = replace_fail(iters_global_2,  maxit); % Replace failures for global with a maxit indicator
    yE2 = replace_fail(iters_sec_rep_2, maxit); % Replace failures for Secant with a maxit indicator
else
    yA2 = iters_near0_2; yB2 = iters_between_2; yC2 = iters_near_r2_2; yD2 = iters_global_2; yE2 = iters_sec_rep_2; % Use original counts if no failures
end

figure(7); clf; hold on; grid on; % Create a new figure for plotting and set up the grid
title('Task 2.3(2): Iterations as rho1,rho2 -> 0 together'); % Set the title of the plot
xlabel('rho1 (log scale), rho2 = 2*rho1'); % Label the x-axis
if hasFail_232 % If there were any failures
    ylabel(sprintf('iterations (failures shown at %d)', maxit+1)); % Label the y-axis indicating failures
else
    ylabel('iterations'); % Label the y-axis normally
end
set(gca,'XScale','log'); % Set the x-axis to logarithmic scale
plot(rho1_list2, yA2, '-o', 'DisplayName','Newton: near0'); % Plot iterations for Newton's method near 0 in the second case
plot(rho1_list2, yB2, '-o', 'DisplayName','Newton: between'); % Plot iterations for Newton's method between roots in the second case
plot(rho1_list2, yC2, '-o', 'DisplayName','Newton: near rho2'); % Plot iterations for Newton's method near r2 in the second case
plot(rho1_list2, yD2, '-o', 'DisplayName','Newton: global'); % Plot iterations for Newton's method global in the second case
plot(rho1_list2, yE2, '-o', 'DisplayName','Secant (rep pair)'); % Plot iterations for Secant method in the second case
legend('Location','northwest'); % Add a legend to the plot
hold off; % Release the hold on the current figure

% ---- distance history (rho2 fixed)
rho2_fixed = 1.0; % Set a fixed value for rho2
rho1_hist  = [1e-2, 1e-6, 1e-10]; % Define a list of rho1 values to test
K_hist = 25; % Set the number of iterations to keep track of

figure(8); clf; hold on; grid on; % Create a new figure for plotting and prepare the grid
title('Task 2.3: Newton distance-to-nearest-root histories (rho2 fixed)'); % Set the title for the plot
xlabel('iteration k'); % Label the x-axis as 'iteration k'
ylabel('distance to nearest root (log scale)'); % Label the y-axis indicating distance to the nearest root
set(gca,'YScale','log'); % Set the y-axis to a logarithmic scale for better visualization

for j = 1:numel(rho1_hist) % Loop through each value of rho1 in the defined list
    r1 = rho1_hist(j); % Get the current rho1 value from the list
    r2 = rho2_fixed; % Use the fixed value of rho2

    fhat  = @(x) x.*(x-r1).*(x-r2); % Define the function fhat based on the current r1 and fixed r2
    fhatp = @(x) 3*x.^2 - 2*(r1+r2)*x + (r1*r2); % Define the derivative of fhat

    x0 = 10*r1;  % Set an initial guess far enough from the root to avoid immediate convergence
    xhist = newton_fixed_steps(fhat, fhatp, x0, K_hist); % Run the fixed-step Newton method and store the history of x values

    roots_now = [0, r1, r2]; % Define the current roots for distance calculation
    dist = zeros(size(xhist)); % Initialize an array to hold the distances to the nearest root
    for k = 1:numel(xhist) % Loop through each value in the history of x
        dist(k) = min(abs(xhist(k) - roots_now)); % Calculate the distance to the nearest root
    end
    dist = max(dist, realmin); % Ensure that distances are not less than the smallest positive number

    plot(0:(numel(dist)-1), dist, '-o', 'DisplayName', sprintf('rho1=%g, x0=10*rho1', r1)); % Plot the distance history with a label for the current rho1
end

legend('Location','northeast'); % Add a legend to the plot in the northeast corner
hold off; % Release the hold on the current figure

fprintf('\nDone. Figures 1–8 generated.\n'); % Print a message indicating that the figures have been generated


