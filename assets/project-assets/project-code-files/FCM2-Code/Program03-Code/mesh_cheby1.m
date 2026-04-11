function x = mesh_cheby1(a, b, n)
% mesh_cheby1
%   Create n+1 Chebyshev points of the first kind on [a, b].
%   Returns the points in a column vector sorted from left to right.
%
% Inputs:
%   a  - left end of the interval
%   b  - right end of the interval
%   n  - degree (number of points = n+1)
%
% Output:
%   x  - column vector of n+1 Chebyshev points on [a, b], increasing order

    x = zeros(n+1, 1); % prepare output vector

    mid   = (a + b) / 2; % midpoint for mapping from [-1,1] to [a,b]
    half  = (b - a) / 2; % half-length of the interval

    for i = 0 : n % go through each index
        t      = cos( (2*i + 1) * pi / (2*(n+1)) ); % chebyshev point on [-1,1]
        x(i+1) = mid + half * t;                     % map to [a,b] and store
    end

    x = sort(x); % ensure points are in increasing order

end
