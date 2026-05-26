function [Q, n_evals] = comp_trap(f, a, b, M)
% Composite trapezoidal rule on [a, b] with M equal subintervals.
% This function uses a one-pass approach: evaluate the integrand at the
% two endpoints and at each interior node once, then combine them to get
% the trapezoidal approximation. Returns the approximation Q and the
% number of function evaluations n_evals (which is M+1).

    H = (b - a) / M;   % compute the width of each subinterval

    % -- evaluate the function at the left endpoint --
    f_a = f(a);        % call f at x = a and store result in f_a

    % -- evaluate the function at the right endpoint --
    f_b = f(b);        % call f at x = b and store result in f_b

    % -- sum the function values at the interior nodes x_1 ... x_{M-1} --
    interior_sum = 0;  % initialize running total for interior node values
    for k = 1 : M-1    % loop over the M-1 interior nodes
        x_k = a + k * H;          % compute the position of the k-th interior node
        interior_sum = interior_sum + f(x_k); % add f(x_k) to the running total
    end

    % -- form the trapezoidal approximation using endpoint and interior sums --
    Q = (H / 2) * (f_a + 2 * interior_sum + f_b); % combine values to get integral estimate

    % -- record the number of function evaluations used --
    n_evals = M + 1;   % endpoints (2) plus M-1 interior nodes gives M+1 evaluations

end