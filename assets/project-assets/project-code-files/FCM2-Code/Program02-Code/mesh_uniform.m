function x = mesh_uniform(a, b, n)
% mesh_uniform
%   Generate n+1 uniformly spaced mesh points on the interval [a, b].
%
%       x_i = a + i*(b-a)/n,  i = 0, 1, ..., n
%
% Inputs:
%   a  - left endpoint of the interval
%   b  - right endpoint of the interval
%   n  - degree of the interpolating polynomial (number of points = n+1)
%
% Outputs:
%   x  - column vector of n+1 uniformly spaced points in [a, b]
%
% Space:  O(n)    Time:  O(n)

    x = zeros(n+1, 1); % Pre-allocate column vector for the mesh points

    h = (b - a) / n; % Compute the uniform step size

    for i = 0 : n % Loop over each mesh index
        x(i+1) = a + i * h; % Compute and store the i-th mesh point (1-indexed)
    end

end
