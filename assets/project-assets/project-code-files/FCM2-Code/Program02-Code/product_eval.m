function p = product_eval(x, alpha, rho)
% product_eval
%   Evaluate a polynomial in product form at each point in x (double precision).
%
%       p_n(x) = alpha * (x - rho(1)) * (x - rho(2)) * ... * (x - rho(n))
%
%   Uses the incremental product algorithm:
%       d_0 = alpha
%       d_i = d_{i-1} * (x - rho(i)),  i = 1, ..., n
%       p_n(x) = d_n
%
%   This algorithm computes p_n(x) to high relative accuracy (Higham 2002),
%   with error |d_n - p_n(x)| <= gamma_{2n+1} * |p_n(x)|,
%   where gamma_k = k*u / (1 - k*u) and u is the unit roundoff.
%   Used as the "exact" denominator when computing condition numbers.
%
% Inputs:
%   x      - vector of query points (double precision)
%   alpha  - leading scalar coefficient of the polynomial
%   rho    - vector of n roots  rho(1), ..., rho(n)
%
% Outputs:
%   p      - vector of polynomial values p_n(x), same size as x
%
% Space:  O(length(x))    Time:  O(n * length(x))

    p = alpha * ones(size(x)); % Initialize running product to alpha at all query points

    for i = 1 : length(rho)        % Loop over each root
        p = p .* (x - rho(i));     % Multiply in the (x - rho_i) factor at all points
    end

end
