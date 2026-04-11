function p = newton_eval_dp(x_nodes, c, x)
% newton_eval_dp
%   Same as newton_eval but works in double precision throughout.
%   We need this for P3 since piecewise_interp and the spline codes all
%   operate in double -- the single-precision casts in newton_eval would
%   silently kill about 8 digits of accuracy.
%
% Inputs:
%   x_nodes - vector of mesh points (double)
%   c       - vector of Newton divided-difference coefficients (double)
%   x       - scalar or vector of query points (double)
%
% Output:
%   p       - interpolated values at x (double)

    n = length(c) - 1;      % polynomial degree

    p = c(n+1) * ones(size(x));   % start from the leading coefficient

    for j = n : -1 : 1            % Horner's method, stepping down
        p = c(j) + (x - x_nodes(j)) .* p;
    end

end