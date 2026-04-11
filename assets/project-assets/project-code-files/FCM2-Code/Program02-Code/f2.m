function y = f2(x, d)
% f2
%   Evaluate f2(x; d) = prod_{i=1}^{d} (x - i) (single precision).
%
%   The product form is available via product_eval with alpha=1
%   and rho = [1, 2, ..., d].  The monomial form is available via
%   integer coefficient recursion:
%       f2(x;1) = (x-1)
%       f2(x;k) = f2(x;k-1) * (x-k)
%
% Inputs:
%   x    - vector of query points (single)
%   d    - degree of the polynomial (number of roots = d)
%
% Outputs:
%   y    - vector of function values prod_{i=1}^{d}(x-i),
%          same size as x (single)
%
% Space:  O(length(x))    Time:  O(d * length(x))

    x = single(x);              % Ensure query points are in single precision
    y = ones(size(x), 'single'); % Initialize running product to 1 at all query points

    for i = 1 : d               % Loop over each integer root
        y = y .* (x - single(i)); % Multiply in the (x - i) factor at all query points
    end

end
