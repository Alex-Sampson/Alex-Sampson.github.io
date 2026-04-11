function val = spline1_eval(s2, breaks, f_vals, x)
% spline1_eval
%   Evaluate the cubic spline s(x) at a single point x, given the second
%   derivatives s2 computed by spline1_coeffs.
%
%   On each subinterval [x_{i-1}, x_i] the local cubic is reconstructed
%   directly from the s_i'' values using the formula derived in Set 13:
%
%     p_i(x) = s2(i-1)/(6*h_i) * (x_i - x)^3
%            + s2(i)  /(6*h_i) * (x - x_{i-1})^3
%            + gamma  * (x - x_{i-1})
%            + gamma_tilde
%
%   where the integration constants gamma and gamma_tilde are recovered
%   from the interpolation conditions at the two endpoints of the interval.
%
% Inputs:
%   s2     - (n+1) x 1 vector of second derivatives at each knot (from spline1_coeffs)
%   breaks - (n+1) x 1 vector of knot locations (same as x_nodes in spline1_coeffs)
%   f_vals - (n+1) x 1 vector of function values at the knots
%   x      - scalar point at which to evaluate s(x)
%
% Output:
%   val    - scalar value s(x)

    n = length(breaks) - 1;   % number of subintervals

    % ---- binary search: find subinterval i such that breaks(i) <= x < breaks(i+1) ----
    % same logic as piecewise_eval -- the rightmost point belongs to the last interval
    lo = 1;
    hi = n;

    while lo < hi
        mid = floor((lo + hi) / 2);
        if x < breaks(mid + 1)
            hi = mid;
        else
            lo = mid + 1;
        end
    end

    i = lo;   % 1-based subinterval index

    % clamp just in case x == breaks(end) lands one past the end
    if i > n
        i = n;
    end

    % ---- grab the local data for subinterval i = [x_{i-1}, x_i] ----
    % note: breaks(i) = x_{i-1} and breaks(i+1) = x_i in 1-based indexing
    x_left  = breaks(i);      % x_{i-1}
    x_right = breaks(i+1);    % x_i
    hi_     = x_right - x_left;   % h_i = x_i - x_{i-1}  (renamed to avoid clash with search var)

    f_left  = f_vals(i);      % f_{i-1}
    f_right = f_vals(i+1);    % f_i

    s2_left  = s2(i);         % s''(x_{i-1})
    s2_right = s2(i+1);       % s''(x_i)

    % ---- recover the two integration constants from Set 13 equations (1)-(2) ----
    % gamma_tilde = f_{i-1} - h_i^2/6 * s''_{i-1}
    gamma_tilde = f_left - (hi_^2 / 6) * s2_left;

    % gamma = (f_i - f_{i-1})/h_i - h_i/6 * (s''_i - s''_{i-1})
    gamma = (f_right - f_left) / hi_ - (hi_ / 6) * (s2_right - s2_left);

    % ---- evaluate the local cubic at x ----
    val = (s2_left  / (6 * hi_)) * (x_right - x)^3 ...   % left cubic term
        + (s2_right / (6 * hi_)) * (x - x_left)^3  ...   % right cubic term
        + gamma       * (x - x_left)               ...   % linear term
        + gamma_tilde;                                     % constant term

end