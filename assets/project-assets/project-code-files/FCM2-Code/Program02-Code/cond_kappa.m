function [kappa_y, kappa_1] = cond_kappa(x_nodes, f_vals, gamma, x, denom)
% cond_kappa
%   Evaluate the interpolation condition numbers kappa(x,n,y) and
%   kappa(x,n,1) at each query point x (double precision).
%
%       kappa(x,n,y) = sum_{i=0}^{n} |l_i(x) * y_i| / |p_n(x)|
%       kappa(x,n,1) = sum_{i=0}^{n} |l_i(x)|
%
%   The numerator of kappa(x,n,y) is computed via Barycentric Form 1
%   with absolute values, using bary1_eval_dp.  The denominator |p_n(x)|
%   is supplied by the caller as denom, which should be the product form
%   evaluation from product_eval when available, or the |p_n(x)| output
%   of bary1_eval_dp otherwise.
%
%   kappa(x,n,1) is the Lebesgue function evaluated at x and does not
%   require a denominator.
%
% Inputs:
%   x_nodes - vector of n+1 distinct mesh points (double)
%   f_vals  - vector of n+1 function values f(x_i) (double)
%   gamma   - vector of n+1 barycentric Form 1 weights (double)
%   x       - vector of query points (double)
%   denom   - vector of "exact" |p_n(x)| values at each query point,
%             from product_eval or bary1_eval_dp (double)
%
% Outputs:
%   kappa_y - vector of kappa(x,n,y) at each query point (double)
%   kappa_1 - vector of kappa(x,n,1) at each query point (double)
%
% Space:  O(n)    Time:  O(n * length(x))

    x_nodes = double(x_nodes); % Ensure nodes are in double precision
    f_vals  = double(f_vals);  % Ensure function values are in double precision
    gamma   = double(gamma);   % Ensure weights are in double precision
    x       = double(x);       % Ensure query points are in double precision
    denom   = double(denom);   % Ensure denominator values are in double precision

    kappa_y = zeros(size(x)); % Pre-allocate kappa(x,n,y) output
    kappa_1 = zeros(size(x)); % Pre-allocate kappa(x,n,1) output

    for k = 1 : length(x)     % Loop over each query point

        diff = x(k) - x_nodes; % Compute differences x - x_i for all nodes

        exact = find(diff == 0); % Check if query point coincides with any node

        if ~isempty(exact)          % If query point is exactly a node
            kappa_y(k) = 1.0;       % Condition number is 1 at an interpolation node
            kappa_1(k) = 1.0;       % Lebesgue function is also 1 at a node
            continue;               % Move to the next query point
        end

        omega = prod(diff);                        % Compute node polynomial omega(x)
        terms = gamma ./ diff;                     % Compute gamma_i / (x - x_i)
        li    = abs(omega) * abs(terms);           % Compute |l_i(x)| for all nodes

        numer_y    = sum(li .* abs(f_vals));       % Absolute numerator of kappa(x,n,y)
        kappa_1(k) = sum(li);                      % kappa(x,n,1) is the Lebesgue function

        if denom(k) == 0                           % Guard against division by zero
            kappa_y(k) = Inf;                      % Condition number is infinite at a root
        else
            kappa_y(k) = numer_y / denom(k);       % Evaluate kappa(x,n,y)
        end

    end

end
