function [p, numer_abs, denom] = bary1_eval_dp(x_nodes, f_vals, gamma, x)
% bary1_eval_dp
%   Evaluate the interpolating polynomial at query points x using
%   Barycentric Form 1 (double precision).
%
%       p_n(x) = omega(x) * sum_{i=0}^{n} gamma_i * f_i / (x - x_i)
%
%   where omega(x) = prod_{i=0}^{n} (x - x_i) is the node polynomial.
%   Also returns the absolute numerator of kappa(x,n,y) and the
%   "exact" denominator |p_n(x)| for use in condition number evaluation.
%
%       numer_abs = sum_{i=0}^{n} |gamma_i * f_i / (x - x_i)| * |omega(x)|
%       denom     = |p_n(x)|  (this output, computed in double precision)
%
% Inputs:
%   x_nodes - vector of n+1 distinct mesh points (double)
%   f_vals  - vector of n+1 function values f(x_i) (double)
%   gamma   - vector of n+1 barycentric Form 1 weights (double)
%   x       - vector of query points at which to evaluate p_n (double)
%
% Outputs:
%   p         - vector of interpolated values, same size as x (double)
%   numer_abs - vector of |numerator of kappa(x,n,y)|, same size as x
%   denom     - vector of |p_n(x)|, same size as x
%
% Space:  O(n)    Time:  O(n * length(x))

    x_nodes = double(x_nodes); % Ensure nodes are in double precision
    f_vals  = double(f_vals);  % Ensure function values are in double precision
    gamma   = double(gamma);   % Ensure weights are in double precision
    x       = double(x);       % Ensure query points are in double precision

    p         = zeros(size(x)); % Pre-allocate interpolated values
    numer_abs = zeros(size(x)); % Pre-allocate absolute numerator of kappa
    denom     = zeros(size(x)); % Pre-allocate denominator |p_n(x)|

    for k = 1 : length(x)      % Loop over each query point

        diff = x(k) - x_nodes; % Compute differences x - x_i for all nodes

        exact = find(diff == 0); % Check if query point coincides with any node

        if ~isempty(exact)                   % If query point is exactly a node
            p(k)         = f_vals(exact(1)); % Return function value directly
            numer_abs(k) = abs(f_vals(exact(1))); % Absolute numerator is just |f_i|
            denom(k)     = abs(f_vals(exact(1))); % Denominator equals |f_i| at a node
            continue;                         % Move to the next query point
        end

        omega = prod(diff);                  % Compute node polynomial omega(x)

        terms     = gamma ./ diff;           % Compute gamma_i / (x - x_i) for all nodes
        p(k)      = omega * sum(f_vals .* terms); % Evaluate Barycentric Form 1

        numer_abs(k) = abs(omega) * sum(abs(f_vals .* terms)); % Absolute numerator of kappa
        denom(k)     = abs(p(k));            % Denominator is |p_n(x)|

    end

end
