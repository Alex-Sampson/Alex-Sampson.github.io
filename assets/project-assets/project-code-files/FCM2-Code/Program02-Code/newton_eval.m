function p = newton_eval(x_nodes, c, x)
% newton_eval
%   Evaluate the Newton form of the interpolating polynomial at query
%   points x using Horner's rule (single precision).
%
%       p_n(x) = c(1) + c(2)*(x-x_1) + c(3)*(x-x_1)*(x-x_2) + ...
%
%   Horner's rule evaluates this as:
%       p_n(x) = c(1) + (x-x_1)*( c(2) + (x-x_2)*( c(3) + ... ))
%
% Inputs:
%   x_nodes - vector of n+1 mesh points in the same order used to
%             compute c (single)
%   c       - vector of n+1 divided difference coefficients from
%             newton_coeffs (single)
%   x       - vector of query points at which to evaluate p_n (single)
%
% Outputs:
%   p       - vector of interpolated values, same size as x (single)
%
% Space:  O(n)    Time:  O(n * length(x))

    x_nodes = single(x_nodes); % Ensure nodes are in single precision
    c       = single(c);       % Ensure coefficients are in single precision
    x       = single(x);       % Ensure query points are in single precision

    n = length(c) - 1;         % Degree of the polynomial

    p = c(n+1) * ones(size(x), 'single'); % Initialize with leading coefficient at all query points

    for j = n : -1 : 1                    % Apply Horner's rule from the highest term downward
        p = c(j) + (x - x_nodes(j)) .* p; % Accumulate the next term at all query points
    end

end
