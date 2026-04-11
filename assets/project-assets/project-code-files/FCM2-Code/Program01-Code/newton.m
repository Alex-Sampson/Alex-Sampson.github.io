function [gamma, f_vals] = bary1_coeffs(x_nodes, n, f_func)
% BARY1_COEFFS  Compute barycentric Form 1 weights (gamma_i) and function values.
%
% Inputs:
%   x_nodes  - vector of n+1 distinct mesh points x(1)..x(n+1)  (1-indexed)
%   n        - degree of polynomial  (number of nodes = n+1)
%   f_func   - function handle, e.g.  @(x) sin(x)
%
% Outputs:
%   gamma    - vector of n+1 barycentric weights  gamma_i = 1 / prod_{j~=i}(x_i - x_j)
%   f_vals   - vector of function values  f(x_i)
%
% Space:  O(n)    Time:  O(n^2)

    gamma  = zeros(n+1, 1);
    f_vals = zeros(n+1, 1);

    for i = 1 : n+1
        f_vals(i) = f_func(x_nodes(i));
        prod_val = 1.0;
        for j = 1 : n+1
            if j ~= i
                prod_val = prod_val * (x_nodes(i) - x_nodes(j));
            end
        end
        gamma(i) = 1.0 / prod_val;
    end
end