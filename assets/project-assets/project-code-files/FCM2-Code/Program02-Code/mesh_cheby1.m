function x = mesh_cheby1(a, b, n)
% mesh_cheby1
%   Generate n+1 Chebyshev points of the first kind on the interval [a, b].
%
%   The Chebyshev points of the first kind are the zeros of the Chebyshev
%   polynomial T_{n+1}(t) on [-1, 1], mapped to [a, b]:
%
%       t_i = cos( (2i+1)*pi / (2*(n+1)) ),  i = 0, 1, ..., n
%       x_i = (a+b)/2 + (b-a)/2 * t_i
%
%   Points are returned in increasing order (left to right).
%
% Inputs:
%   a  - left endpoint of the interval
%   b  - right endpoint of the interval
%   n  - degree of the interpolating polynomial (number of points = n+1)
%
% Outputs:
%   x  - column vector of n+1 Chebyshev first kind points in [a, b],
%        sorted in increasing order
%
% Space:  O(n)    Time:  O(n)

    x = zeros(n+1, 1); % Pre-allocate column vector for the mesh points

    mid   = (a + b) / 2; % Compute the midpoint of the interval for the affine map
    half  = (b - a) / 2; % Compute the half-width of the interval for the affine map

    for i = 0 : n % Loop over each Chebyshev index
        t      = cos( (2*i + 1) * pi / (2*(n+1)) ); % Compute the i-th Chebyshev point on [-1, 1]
        x(i+1) = mid + half * t;                     % Map the point to [a, b] and store (1-indexed)
    end

    x = sort(x); % Sort points into increasing order for consistent downstream use

end
