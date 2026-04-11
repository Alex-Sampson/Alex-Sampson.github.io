function val = spline2_eval(alpha, breaks, x)
% spline2_eval
%   Evaluate the B-spline cubic spline s(x) at a single point x.
%
%   s(x) = sum_{i=-1}^{n+1} alpha_i * B_i(x)
%
%   On any subinterval [x_{j}, x_{j+1}] exactly four B-splines are nonzero:
%   B_{j-1}, B_j, B_{j+1}, B_{j+2} (from the support property in Set 14).
%   We evaluate each of those four using the explicit piecewise definition
%   from Definition 3.1 and accumulate the sum.
%
% Inputs:
%   alpha  - (n+3) x 1 vector of B-spline coefficients from spline2_coeffs
%            stored as [alpha_{-1}; alpha_0; ...; alpha_{n+1}]
%   breaks - (n+1) x 1 vector of uniform knot locations
%   x      - scalar query point
%
% Output:
%   val    - scalar value s(x)

    n = length(breaks) - 1;
    h = breaks(2) - breaks(1);   % uniform spacing

    % ---- binary search: find j such that breaks(j) <= x < breaks(j+1) ----
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

    j = lo - 1;   % convert to 0-based subinterval index [x_j, x_{j+1}]

    % clamp for the rightmost endpoint
    if j >= n
        j = n - 1;
    end

    % ---- evaluate the four active B-splines on [x_j, x_{j+1}] ----
    % the active ones are B_{j-1}, B_j, B_{j+1}, B_{j+2}
    % each B_i is defined by Definition 3.1 in Set 14 -- we evaluate
    % by computing the local variable t = x - x_i and checking which
    % piece of the piecewise definition applies.
    %
    % storage convention: alpha_i is at position i+2 (so alpha_{-1} is at 1)

    val = 0;

    for k = j-1 : j+2   % loop over the four active B-spline indices (0-based)

        % reference center of B_k is x_k = breaks(k+1) in 1-based storage
        % but k can be -1 or n+1 which are outside breaks -- that's fine,
        % we just compute the center by extrapolation from the uniform mesh
        x_k = breaks(1) + k * h;   % center of B_k (uniform mesh so this is exact)

        % local distance from the center, normalized to units of h
        % the five-point support of B_k is [x_k - 2h, x_k + 2h]
        t = x - x_k;   % raw distance

        % evaluate B_k(x) using Definition 3.1 piece by piece
        % (the definition uses the variable (x - x_{i-2}) etc; here t = x - x_k)
        if t <= -2*h || t >= 2*h
            Bval = 0;   % outside support

        elseif t <= -h
            % piece 1: x in [x_{k-2}, x_{k-1}], i.e. t in [-2h, -h]
            % B_k(x) = (1/h^3) * (x - x_{k-2})^3 = (1/h^3) * (t + 2h)^3
            Bval = (t + 2*h)^3 / h^3;

        elseif t <= 0
            % piece 2: x in [x_{k-1}, x_k], i.e. t in [-h, 0]
            % B_k(x) = (1/h^3) * (h^3 + 3h^2*(t+h) + 3h*(t+h)^2 - 3*(t+h)^3)
            u = t + h;   % u = x - x_{k-1}, ranges in [0, h]
            Bval = (h^3 + 3*h^2*u + 3*h*u^2 - 3*u^3) / h^3;

        elseif t <= h
            % piece 3: x in [x_k, x_{k+1}], i.e. t in [0, h]
            % B_k(x) = (1/h^3) * (h^3 + 3h^2*(h-t) + 3h*(h-t)^2 - 3*(h-t)^3)
            u = h - t;   % u = x_{k+1} - x, ranges in [0, h]
            Bval = (h^3 + 3*h^2*u + 3*h*u^2 - 3*u^3) / h^3;

        else
            % piece 4: x in [x_{k+1}, x_{k+2}], i.e. t in [h, 2h]
            % B_k(x) = (1/h^3) * (x_{k+2} - x)^3 = (1/h^3) * (2h - t)^3
            Bval = (2*h - t)^3 / h^3;
        end

        % accumulate: alpha_k is at storage position k+2
        val = val + alpha(k + 2) * Bval;

    end

end