function c = newton_coeffs(x_nodes, f_vals)
% newton_coeffs
%   Compute the divided difference coefficients for the Newton form
%   of the interpolating polynomial (single precision).
%
%       p_n(x) = c(1) + c(2)*(x-x_1) + c(3)*(x-x_1)*(x-x_2) + ...
%
%   Coefficients are computed in-place using the standard recursion:
%
%       c(i) = [x_1, ..., x_i] f  (i-th divided difference)
%
% Inputs:
%   x_nodes - vector of n+1 distinct mesh points, in desired order (single)
%   f_vals  - vector of n+1 function values f(x_i) (single)
%
% Outputs:
%   c       - vector of n+1 divided difference coefficients (single)
%
% Space:  O(n)    Time:  O(n^2)

    x_nodes = single(x_nodes); % Ensure nodes are in single precision
    c       = single(f_vals);  % Copy f_vals into c; will be overwritten in-place

    n = length(c) - 1;         % Degree of the polynomial

    for j = 1 : n              % Loop over each level of the divided difference table
        for i = n : -1 : j    % Update entries in reverse to allow in-place computation
            c(i+1) = (c(i+1) - c(i)) / (x_nodes(i+1) - x_nodes(i-j+1)); % Compute divided difference
        end
    end

end
