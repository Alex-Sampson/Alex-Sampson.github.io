% task2driver_p3.m
%
% Task 2: Real data application.
%
% Given 8 discrete data points (t_i, y_i) on a nonuniform mesh, fit
% interpolants and use them to recover three related functions:
%
%   y(t)  -- the interpolant itself
%   f(t)  -- recovered via f(t) = y(t) + t * y'(t)
%   D(t)  -- recovered via D(t) = exp(-t) * y(t)
%
% Part 1: Natural cubic spline s(t) through the data (nonuniform mesh).
%   s'(t) is obtained analytically by differentiating the local cubic
%   reconstruction formula from Set 13 -- this is exact calculus on the
%   spline, not numerical differentiation.
%
% Part 2: Piecewise polynomial g_d(t) investigation.
%   Standard piecewise Lagrange (d=1,2,3) is only C^0. Its derivative
%   is discontinuous at the breakpoints, so f(t) = y + t*y' cannot be
%   estimated reliably -- the derivative jumps at every knot. y(t) and
%   D(t) are fine since they require no derivative.
%
%   Piecewise cubic Hermite is C^1 and avoids this problem. We prescribe
%   endpoint derivatives using s'(t_i) from the spline -- the only
%   consistent source of derivative estimates available from this data.
%
%   All piecewise methods use the exact data knot locations as breakpoints
%   via the updated piecewise_interp (9th argument breaks_in = t_data).
%
% Outputs:
%   Table 1:  Tabulated y(t), f(t), D(t) from spline, t=0.5 to 20, dt=0.5
%   Table 2:  Tabulated y(t), f(t), D(t) from PW Hermite
%   Table 3:  Max |spline - method| for y and D across all methods,
%             and max |spline - Hermite| for f
%   Figure 1: y(t) -- all methods with data points marked
%   Figure 2: f(t) -- spline and PW Hermite only (Lagrange excluded:
%             its discontinuous derivative makes f(t) unreliable)
%   Figure 3: D(t) -- all methods (no derivative needed)
%   Figure 4: s'(t) from spline -- motivates why C^1 matters for f(t)

clear; clc; close all;

% ---- color palette ----
clr_spline = [0.122, 0.471, 0.706];   % steel blue   -- cubic spline
clr_herm   = [0.890, 0.102, 0.110];   % crimson      -- PW Hermite
clr_pw1    = [0.200, 0.627, 0.173];   % forest green -- PW d=1
clr_pw2    = [0.596, 0.306, 0.639];   % purple       -- PW d=2
clr_pw3    = [1.000, 0.498, 0.000];   % orange       -- PW d=3
clr_data   = [0.000, 0.000, 0.000];   % black        -- data points

% ---- the 8 given data points ----
t_data = [0.5;  1.0;  2.0;  4.0;  5.0;  10.0; 15.0; 20.0];
y_data = [0.0552; 0.0600; 0.0682; 0.0801; 0.0843; 0.0931; 0.0912; 0.0857];

n_data = length(t_data);   % 8 points, 7 subintervals
M      = n_data - 1;       % number of subintervals = 7

% ---- tabulation grid: t = 0.5 to 20.0, dt = 0.5 (40 points) ----
dt    = 0.5;
t_tab = (0.5 : dt : 20.0)';
n_tab = length(t_tab);

fprintf('\n============================================================\n');
fprintf('TASK 2: Real Data Application\n');
fprintf('%d data points, tabulating %d points from t=%.1f to t=%.1f\n', ...
    n_data, n_tab, t_tab(1), t_tab(end));
fprintf('============================================================\n\n');

% ============================================================
% PART 1: Natural cubic spline
% ============================================================

fprintf('------------------------------------------------------------\n');
fprintf('PART 1: Natural cubic spline s(t)\n');
fprintf('------------------------------------------------------------\n\n');

% build the spline -- spline1_coeffs supports nonuniform meshes natively
[s2, br, fv] = spline1_coeffs(t_data, y_data, 'natural', []);

% evaluate y(t) = s(t) on the tabulation grid
y_spl = zeros(n_tab, 1);
for k = 1 : n_tab
    y_spl(k) = spline1_eval(s2, br, fv, t_tab(k));
end

% evaluate s'(t) analytically -- differentiate the Set 13 local cubic.
% on subinterval [t_{i-1}, t_i] with h_i = t_i - t_{i-1}:
%
%   p_i(t) = s2(i-1)/(6h)*(t_i-t)^3 + s2(i)/(6h)*(t-t_{i-1})^3
%           + gamma*(t-t_{i-1}) + gamma_tilde
%
%   p_i'(t) = -s2(i-1)/(2h)*(t_i-t)^2 + s2(i)/(2h)*(t-t_{i-1})^2 + gamma
%
% where gamma = (f_i - f_{i-1})/h - h/6*(s2(i) - s2(i-1))
dy_spl = zeros(n_tab, 1);
for k = 1 : n_tab
    dy_spl(k) = spline_deriv(t_tab(k), s2, br, fv, n_data);
end

% recover f(t) and D(t) from the spline
f_spl = y_spl + t_tab .* dy_spl;    % f(t) = y(t) + t * s'(t)
D_spl = exp(-t_tab) .* y_spl;       % D(t) = exp(-t) * y(t)

% print Table 1
fprintf('Table 1: Spline estimates of y(t), f(t), D(t)\n');
fprintf('%-8s %-14s %-14s %-14s\n', 't', 'y(t)', 'f(t)', 'D(t)');
fprintf('%s\n', repmat('-', 1, 52));
for k = 1 : n_tab
    fprintf('%-8.1f %-14.6f %-14.6f %-14.6f\n', ...
        t_tab(k), y_spl(k), f_spl(k), D_spl(k));
end
fprintf('\n');

% ============================================================
% PART 2: Piecewise polynomial investigation
% ============================================================

fprintf('------------------------------------------------------------\n');
fprintf('PART 2: Piecewise polynomial g_d(t)\n');
fprintf('------------------------------------------------------------\n\n');

fprintf('Standard piecewise Lagrange (d=1,2,3) is only C^0.\n');
fprintf('Its derivative is discontinuous at the 6 interior breakpoints.\n');
fprintf('f(t) = y(t) + t*y''(t) CANNOT be estimated reliably from these\n');
fprintf('methods -- the derivative jumps at each knot.\n');
fprintf('y(t) and D(t) are estimated since they need no derivative.\n\n');
fprintf('PW Hermite is C^1 and gives a smooth derivative everywhere.\n');
fprintf('Endpoint derivatives are prescribed using s''(t_i) from the spline.\n\n');

% use the spline as a continuous function handle -- this lets piecewise_interp
% evaluate f at any point including interior nodes within each subinterval.
% at the data knots it returns the exact data values by spline interpolation.
f_handle = @(t) spline1_eval(s2, br, fv, t);

% compute spline derivatives at each data knot for the Hermite BC
dy_at_knots = zeros(n_data, 1);
for i = 1 : n_data
    dy_at_knots(i) = spline_deriv(t_data(i), s2, br, fv, n_data);
end

% derivative handle -- evaluated at the knots only (endpoints of each subinterval)
% piecewise_interp calls f_prime at a_s and b_s which are always data knots,
% so linear interpolation between knots is exact here
df_handle = @(t) interp1(t_data, dy_at_knots, t, 'linear', 'extrap');

% build all piecewise interpolants on the exact data knot mesh.
% passes t_data as breaks_in (9th argument) so the subinterval boundaries
% match the data knot locations exactly -- not a uniform approximation.
[coef_d1, bk_d1, nd_d1] = piecewise_interp(f_handle, [],        t_data(1), t_data(end), M, 1, 'uniform', false, t_data);
[coef_d2, bk_d2, nd_d2] = piecewise_interp(f_handle, [],        t_data(1), t_data(end), M, 2, 'uniform', false, t_data);
[coef_d3, bk_d3, nd_d3] = piecewise_interp(f_handle, [],        t_data(1), t_data(end), M, 3, 'uniform', false, t_data);
[coef_ph, bk_ph, nd_ph] = piecewise_interp(f_handle, df_handle, t_data(1), t_data(end), M, 3, 'uniform', true,  t_data);

% evaluate y(t) on the tabulation grid
y_d1   = zeros(n_tab,1); y_d2   = zeros(n_tab,1);
y_d3   = zeros(n_tab,1); y_herm = zeros(n_tab,1);

for k = 1 : n_tab
    y_d1(k)   = piecewise_eval(coef_d1, bk_d1, nd_d1, t_tab(k), false);
    y_d2(k)   = piecewise_eval(coef_d2, bk_d2, nd_d2, t_tab(k), false);
    y_d3(k)   = piecewise_eval(coef_d3, bk_d3, nd_d3, t_tab(k), false);
    y_herm(k) = piecewise_eval(coef_ph, bk_ph, nd_ph, t_tab(k), true);
end

% for PW Hermite, recover y'(t) analytically from the local Newton polynomial.
% the Hermite cubic on each subinterval is stored in Newton form with repeated
% nodes [a_s, a_s, b_s] -- we differentiate it by using the fact that the
% derivative of the Newton polynomial with nodes [x0, x1, x2] evaluated
% at t is: c2 + c3*((t-x0) + (t-x1)) using the product rule on Horner form.
% this is exact calculus on the stored polynomial, not an approximation.
dy_herm = zeros(n_tab, 1);
for k = 1 : n_tab
    dy_herm(k) = hermite_deriv(t_tab(k), coef_ph, bk_ph, nd_ph, M);
end

% recover f(t) and D(t)
f_herm = y_herm + t_tab .* dy_herm;

D_d1   = exp(-t_tab) .* y_d1;
D_d2   = exp(-t_tab) .* y_d2;
D_d3   = exp(-t_tab) .* y_d3;
D_herm = exp(-t_tab) .* y_herm;

% print Table 2: PW Hermite
fprintf('Table 2: PW Hermite estimates of y(t), f(t), D(t)\n');
fprintf('%-8s %-14s %-14s %-14s\n', 't', 'y(t)', 'f(t)', 'D(t)');
fprintf('%s\n', repmat('-', 1, 52));
for k = 1 : n_tab
    fprintf('%-8.1f %-14.6f %-14.6f %-14.6f\n', ...
        t_tab(k), y_herm(k), f_herm(k), D_herm(k));
end
fprintf('\n');

% print Table 3
fprintf('Table 3: Max |spline - method| over %d tabulation points\n', n_tab);
fprintf('%-18s %-14s %-14s\n', 'Method', 'max |Dy|', 'max |DD|');
fprintf('%s\n', repmat('-', 1, 48));
pw_names = {'PW d=1', 'PW d=2', 'PW d=3', 'PW Hermite'};
y_pw_all = {y_d1, y_d2, y_d3, y_herm};
D_pw_all = {D_d1, D_d2, D_d3, D_herm};
for m = 1 : 4
    fprintf('%-18s %-14.4e %-14.4e\n', pw_names{m}, ...
        max(abs(y_pw_all{m} - y_spl)), ...
        max(abs(D_pw_all{m} - D_spl)));
end
fprintf('\n');
fprintf('f(t) comparison (C^1 methods only):\n');
fprintf('%-18s %-14s\n', 'Method', 'max |Df|');
fprintf('%s\n', repmat('-', 1, 34));
fprintf('%-18s %-14.4e\n', 'PW Hermite', max(abs(f_herm - f_spl)));
fprintf('\n');

% ============================================================
% FIGURES
% ============================================================

% ---- Figure 1: y(t) -- all methods with data points ----
% all methods pass through the data knots by construction.
% differences appear between the knots, particularly in the wide
% gaps [5,10] and [10,15] of this nonuniform mesh.
figure(1); clf; hold on; grid on; box on;
plot(t_tab, y_spl,  '-',  'Color', clr_spline, 'LineWidth', 2.0, 'DisplayName', 'Spline (natural)');
plot(t_tab, y_herm, '--', 'Color', clr_herm,   'LineWidth', 1.8, 'DisplayName', 'PW Hermite');
plot(t_tab, y_d3,   ':',  'Color', clr_pw3,    'LineWidth', 1.8, 'DisplayName', 'PW d=3');
plot(t_tab, y_d2,   ':',  'Color', clr_pw2,    'LineWidth', 1.8, 'DisplayName', 'PW d=2');
plot(t_tab, y_d1,   ':',  'Color', clr_pw1,    'LineWidth', 1.8, 'DisplayName', 'PW d=1');
plot(t_data, y_data, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'DisplayName', 'Data (t_i, y_i)');
xlabel('t'); ylabel('y(t)');
title('Figure 1: Estimates of y(t) from all methods');
legend('Location', 'best'); hold off;

% ---- Figure 2: f(t) = y(t) + t*y'(t) -- C^1 methods only ----
% standard Lagrange is excluded because its derivative is discontinuous
% at breakpoints, making f(t) unreliable. spline and Hermite give
% smooth, comparable estimates.
figure(2); clf; hold on; grid on; box on;
plot(t_tab, f_spl,  '-',  'Color', clr_spline, 'LineWidth', 2.0, 'DisplayName', 'Spline (natural)');
plot(t_tab, f_herm, '--', 'Color', clr_herm,   'LineWidth', 1.8, 'DisplayName', 'PW Hermite');
plot(t_data, y_data + t_data.*dy_at_knots, 'ko', 'MarkerSize', 6, ...
    'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
xlabel('t'); ylabel('f(t) = y(t) + t \cdot y''(t)');
title('Figure 2: Estimates of f(t) -- C^1 methods only (PW Lagrange excluded)');
legend('Location', 'best'); hold off;

% ---- Figure 3: D(t) = exp(-t)*y(t) -- all methods ----
% no derivative needed so all methods give smooth curves.
% differences reflect only the y(t) interpolation quality.
figure(3); clf; hold on; grid on; box on;
plot(t_tab, D_spl,  '-',  'Color', clr_spline, 'LineWidth', 2.0, 'DisplayName', 'Spline (natural)');
plot(t_tab, D_herm, '--', 'Color', clr_herm,   'LineWidth', 1.8, 'DisplayName', 'PW Hermite');
plot(t_tab, D_d3,   ':',  'Color', clr_pw3,    'LineWidth', 1.8, 'DisplayName', 'PW d=3');
plot(t_tab, D_d2,   ':',  'Color', clr_pw2,    'LineWidth', 1.8, 'DisplayName', 'PW d=2');
plot(t_tab, D_d1,   ':',  'Color', clr_pw1,    'LineWidth', 1.8, 'DisplayName', 'PW d=1');
plot(t_data, exp(-t_data).*y_data, 'ko', 'MarkerSize', 6, ...
    'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
xlabel('t'); ylabel('D(t) = e^{-t} y(t)');
title('Figure 3: Estimates of D(t) -- all methods (no derivative needed)');
legend('Location', 'best'); hold off;

% ---- Figure 4: spline derivative s'(t) ----
% shows the smooth continuous derivative that the spline provides and
% the discrete knot values used to seed the Hermite interpolant.
% makes concrete why C^1 continuity is essential for computing f(t).
figure(4); clf; hold on; grid on; box on;
plot(t_tab, dy_spl, '-', 'Color', clr_spline, 'LineWidth', 2.0, 'DisplayName', 's''(t) -- spline');
plot(t_tab, dy_herm, '--', 'Color', clr_herm, 'LineWidth', 1.8, 'DisplayName', 'y''(t) -- PW Hermite');
plot(t_data, dy_at_knots, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', clr_herm, ...
    'DisplayName', 's''(t_i) used as Hermite BCs');
xline(t_data(2:end-1)', '--', 'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');
xlabel('t'); ylabel('y''(t)');
title('Figure 4: Derivative estimates -- spline and PW Hermite (grey lines mark breakpoints)');
legend('Location', 'best'); hold off;

fprintf('DONE. Produced Tables 1-3 and Figures 1-4.\n');

% ============================================================
% LOCAL HELPER FUNCTIONS
% ============================================================

function dy = spline_deriv(t, s2, br, fv, n_data)
% spline_deriv
%   Analytically evaluate s'(t) by differentiating the Set 13 local cubic.
%   This is exact calculus on the spline -- not a numerical approximation.

    % binary search for the correct subinterval
    lo = 1;  hi = n_data - 1;
    while lo < hi
        mid = floor((lo + hi) / 2);
        if t < br(mid + 1)
            hi = mid;
        else
            lo = mid + 1;
        end
    end
    i = lo;
    if i > n_data - 1,  i = n_data - 1;  end

    t_left   = br(i);       t_right  = br(i+1);
    h_i      = t_right - t_left;
    f_left   = fv(i);       f_right  = fv(i+1);
    s2_left  = s2(i);       s2_right = s2(i+1);

    gamma = (f_right - f_left) / h_i - (h_i / 6) * (s2_right - s2_left);

    dy = -(s2_left  / (2*h_i)) * (t_right - t)^2 ...
         +(s2_right / (2*h_i)) * (t - t_left)^2  ...
         + gamma;
end

function dy = hermite_deriv(t, coeffs, breaks, node_coords, M)
% hermite_deriv
%   Analytically evaluate the derivative of the piecewise cubic Hermite
%   polynomial at a scalar point t.
%
%   The Hermite cubic on each subinterval is stored in Newton form with
%   nodes [a_s, a_s, b_s] and coefficients [c1, c2, c3, c4].
%   The Newton polynomial is:
%     p(t) = c1 + c2*(t-a) + c3*(t-a)^2 + c4*(t-a)^2*(t-b)
%   Differentiating:
%     p'(t) = c2 + 2*c3*(t-a) + c4*(2*(t-a)*(t-b) + (t-a)^2)

    % binary search for the correct subinterval
    lo = 1;  hi = M;
    while lo < hi
        mid = floor((lo + hi) / 2);
        if t < breaks(mid + 1)
            hi = mid;
        else
            lo = mid + 1;
        end
    end
    s = lo;
    if s > M,  s = M;  end

    % local Newton coefficients and repeated nodes for this subinterval
    c1 = coeffs(1, s);
    c2 = coeffs(2, s);
    c3 = coeffs(3, s);
    c4 = coeffs(4, s);

    % nodes stored as [a_s, a_s, b_s, b_s] -- first three used in Horner
    a_s = node_coords(1, s);   % = a_s (repeated)
    b_s = node_coords(3, s);   % = b_s

    % derivative of the Newton form p(t) = c1 + c2*(t-a) + c3*(t-a)^2 + c4*(t-a)^2*(t-b)
    ta = t - a_s;
    tb = t - b_s;

    dy = c2 + 2*c3*ta + c4*(2*ta*tb + ta^2);
end