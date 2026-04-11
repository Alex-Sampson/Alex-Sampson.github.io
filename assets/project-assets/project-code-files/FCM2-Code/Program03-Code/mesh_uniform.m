function x = mesh_uniform(a, b, n)
% mesh_uniform
%   Create n+1 evenly spaced points from a to b.
%
% Inputs:
%   a  - left end of the interval
%   b  - right end of the interval
%   n  - degree of interpolating polynomial (number of points = n+1)
%
% Output:
%   x  - column vector with n+1 evenly spaced points in [a, b]

    x = zeros(n+1, 1); % make column vector to hold the points

    h = (b - a) / n; % step size between consecutive points

    for i = 0 : n % go through each index
        x(i+1) = a + i * h; % compute the i-th point (MATLAB is 1-based)
    end

end
