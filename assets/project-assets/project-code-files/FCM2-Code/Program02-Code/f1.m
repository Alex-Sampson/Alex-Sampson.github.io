function y = f1(x, rho, d)
% f1
%   Evaluate f1(x) = (x - rho)^d (single precision).
%
%   The monomial coefficients are available via the binomial expansion
%   and the product form is available via product_eval with alpha=1
%   and rho repeated d times.
%
% Inputs:
%   x    - vector of query points (single)
%   rho  - root of the polynomial
%   d    - degree of the polynomial
%
% Outputs:
%   y    - vector of function values (x - rho)^d, same size as x (single)
%
% Space:  O(length(x))    Time:  O(length(x))

    y = single((single(x) - single(rho)) .^ d); % Evaluate (x - rho)^d in single precision

end
