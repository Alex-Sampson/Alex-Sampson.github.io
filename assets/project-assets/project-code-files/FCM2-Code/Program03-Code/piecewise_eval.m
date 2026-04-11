function val = piecewise_eval(coeffs, breaks, node_coords, x, hermite_flag)
% piecewise_eval
%   Evaluate the piecewise polynomial g_d(x) at a single point x.
%   All the hard work (building the coefficients) already happened in
%   piecewise_interp -- this just finds the right subinterval and calls
%   newton_eval on it.
%
% Inputs:
%   coeffs       - (d+1) x M matrix of Newton coefficients from piecewise_interp
%   breaks       - (M+1) x 1 vector of subinterval breakpoints
%   node_coords  - (d+1) x M matrix of node locations from piecewise_interp
%   x            - scalar point where we want to evaluate g_d(x)
%   hermite_flag - true if the coefficients came from a Hermite build
%                  (changes which node sequence gets passed to newton_eval)
%
% Output:
%   val - the scalar value g_d(x)

    M = length(breaks) - 1;   % number of subintervals

    % ---- binary search: find which subinterval contains x ----
    % we want the index s such that breaks(s) <= x < breaks(s+1).
    % the rightmost point x == b belongs to the last subinterval.

    lo = 1;
    hi = M;

    while lo < hi
        mid = floor((lo + hi) / 2);    % midpoint of the search range
        if x < breaks(mid + 1)
            hi = mid;                  % x is in the left half
        else
            lo = mid + 1;              % x is in the right half
        end
    end

    s = lo;   % s is now the 1-based index of the correct subinterval

    % clamp to [1, M] just in case x lands exactly on b due to floating point
    if s > M
        s = M;
    end

    % ---- grab the local nodes and coefficients for subinterval s ----
    c       = coeffs(:, s);        % Newton coefficients for this subinterval
    x_nodes = node_coords(:, s);   % node locations for this subinterval

    % for Hermite the Newton evaluation nodes are [a_s, a_s, b_s] (first 3 entries),
    % not all 4 stored coordinates -- degree 3 polynomial needs 3 subtracted nodes
    if hermite_flag
        x_nodes = x_nodes(1:3);
    end

    % ---- evaluate the local Newton polynomial at x ----
    % newton_eval expects column vectors; x is a scalar here
    val = newton_eval_dp(x_nodes, c, x);

end