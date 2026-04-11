function x = mesh_cheby2(a, b, n)
% mesh_cheby2
%   Create n+1 Chebyshev points of the second kind on [a, b].
%   These are the cosine-spaced points that include both interval ends.
%   The function returns a column vector of points sorted from a to b.
%
% Inputs:
%   a  - left endpoint
%   b  - right endpoint
%   n  - polynomial degree (number of points = n+1)
%
% Output:
%   x  - column vector of n+1 Chebyshev second-kind points on [a, b]
%
% Memory/time: linear in n.

    x = zeros(n+1, 1); % prepare output vector

    mid  = (a + b) / 2; % midpoint of the interval
    half = (b - a) / 2; % half the interval length

    for i = 0 : n % go through each index 0..n
        t      = cos( i * pi / n );  % cosine-spaced point on [-1,1]
        x(i+1) = mid + half * t;     % map to [a,b] and store (MATLAB is 1-based)
    end

    x = sort(x); % ensure points are in increasing order

end
