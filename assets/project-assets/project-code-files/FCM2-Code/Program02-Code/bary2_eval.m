function p = bary2_eval(x_nodes, f_vals, beta, x)
% bary2_eval
%   Evaluate the interpolating polynomial at query points x using
%   Barycentric Form 2 (single precision).
%
%       p_n(x) = [ sum_{i=0}^{n} beta_i * f_i / (x - x_i) ]
%                / [ sum_{i=0}^{n} beta_i / (x - x_i) ]
%
%   The weights beta are supplied by bary2_coeffs, which provides O(n)
%   closed-form weights for uniform, cheby1 and cheby2 meshes
%   (Berrut & Trefethen 2004, eqs. 5.1, 5.3, 5.4), or the O(n^2)
%   product formula for general meshes.
%
%   Any common nonzero factor in all beta_i cancels between numerator
%   and denominator, so scaled weights produce the same result.
%
%   If x coincides exactly with a node x_i, the function value f_i
%   is returned directly (Berrut & Trefethen 2004, Section 7).
%
% Inputs:
%   x_nodes - vector of n+1 distinct mesh points (single)
%   f_vals  - vector of n+1 function values f(x_i) (single)
%   beta    - vector of n+1 barycentric Form 2 weights from bary2_coeffs
%             (single)
%   x       - vector of query points at which to evaluate p_n (single)
%
% Outputs:
%   p       - vector of interpolated values, same size as x (single)
%
% Space:  O(n)    Time:  O(n * length(x))

    x_nodes = single(x_nodes); % Ensure nodes are in single precision
    f_vals  = single(f_vals);  % Ensure function values are in single precision
    beta    = single(beta);    % Ensure weights are in single precision
    x       = single(x);       % Ensure query points are in single precision

    p = zeros(size(x), 'single'); % Pre-allocate output in single precision

    for k = 1 : length(x)         % Loop over each query point

        diff = x(k) - x_nodes;    % Compute differences x - x_i for all nodes

        exact = find(diff == 0);   % Check if query point coincides with any node

        if ~isempty(exact)                  % If query point is exactly a node
            p(k) = f_vals(exact(1));        % Return the function value directly
            continue;                        % Move to the next query point
        end

        w       = beta ./ diff;    % Compute beta_i / (x - x_i) for all nodes
        numer   = sum(w .* f_vals); % Compute the numerator sum
        denom   = sum(w);           % Compute the denominator sum

        p(k) = numer / denom;      % Evaluate Barycentric Form 2 at this query point

    end

end
