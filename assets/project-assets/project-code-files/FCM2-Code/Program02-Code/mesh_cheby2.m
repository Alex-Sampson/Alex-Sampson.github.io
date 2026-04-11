function x = mesh_cheby2(a, b, n)
% mesh_cheby2
%   Generate n+1 Chebyshev points of the second kind on the interval [a, b].
%
%   The Chebyshev points of the second kind are the extrema of the Chebyshev
%   polynomial T_n(t) on [-1, 1], mapped to [a, b]:
%
%       t_i = cos( i*pi / n ),  i = 0, 1, ..., n
%       x_i = (a+b)/2 + (b-a)/2 * t_i
%
%   Note: the endpoints t=+1 and t=-1 are included (i=0 and i=n),
%   which distinguishes this mesh from the first kind.
%   Points are returned in increasing order (left to right).
%
% Inputs:
%   a  - left endpoint of the interval
%   b  - right endpoint of the interval
%   n  - degree of the interpolating polynomial (number of points = n+1)
%
% Outputs:
%   x  - column vector of n+1 Chebyshev second kind points in [a, b],
%        sorted in increasing order
%
% Space:  O(n)    Time:  O(n)

    x = zeros(n+1, 1); % Pre-allocate column vector for the mesh points

    mid  = (a + b) / 2; % Compute the midpoint of the interval for the affine map
    half = (b - a) / 2; % Compute the half-width of the interval for the affine map

    for i = 0 : n % Loop over each Chebyshev index
        t      = cos( i * pi / n );  % Compute the i-th extremum of T_n on [-1, 1]
        x(i+1) = mid + half * t;     % Map the point to [a, b] and store (1-indexed)
    end

    x = sort(x); % Sort points into increasing order for consistent downstream use

end
