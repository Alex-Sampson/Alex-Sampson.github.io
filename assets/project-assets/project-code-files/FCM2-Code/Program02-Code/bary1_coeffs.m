function [gamma, f_vals] = bary1_coeffs(x_nodes, n, f_func)
% bary1_coeffs
%   Compute barycentric Form 1 weights and function values (double precision).
%
%       gamma_i = 1 / prod_{j ~= i} (x_i - x_j),  i = 0, 1, ..., n
%
% Inputs:
%   x_nodes - vector of n+1 distinct mesh points (double)
%   n       - degree of the polynomial (number of nodes = n+1)
%   f_func  - function handle for the interpolated function, e.g. @(x) sin(x)
%
% Outputs:
%   gamma   - vector of n+1 barycentric Form 1 weights (double)
%   f_vals  - vector of n+1 function values f(x_i) (double)
%
% Space:  O(n)    Time:  O(n^2)

    gamma  = zeros(n+1, 1); % Pre-allocate barycentric weight vector
    f_vals = zeros(n+1, 1); % Pre-allocate function value vector

    for i = 1 : n+1                      % Loop over each node
        f_vals(i) = f_func(x_nodes(i));  % Evaluate function at the i-th node
        prod_val  = 1.0;                  % Initialize running product for this node

        for j = 1 : n+1          % Loop over all other nodes
            if j ~= i             % Skip the i-th node itself
                prod_val = prod_val * (x_nodes(i) - x_nodes(j)); % Accumulate (x_i - x_j)
            end
        end

        gamma(i) = 1.0 / prod_val; % Weight is the reciprocal of the node product
    end

end
