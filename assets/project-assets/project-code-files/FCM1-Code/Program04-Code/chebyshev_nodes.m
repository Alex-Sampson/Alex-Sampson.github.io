function x = chebyshev_nodes(d)
%    
% CHEBYSHEV_NODES  Roots of the Chebyshev polynomial T_{d+1}(x).
%
% For a nonnegative integer d, the Chebyshev nodes are the roots of
% T_{d+1}(x) and are given by
%
%     x_i = cos( (2*i + 1) * pi / (2*(d+1)) ),   i = 0,1,...,d.
%
% These d+1 points lie in [-1,1].  They are often used as interpolation
% nodes or evaluation points with the Chebyshev basis.
%
% INPUT:
%   d : nonnegative integer (degree)
%
% OUTPUT:
%   x : (d+1)-by-1 vector of Chebyshev nodes in [-1,1]
%
% ----------------------------------------------------------------------------------------------------

    
    if d < 0 || floor(d) ~= d
        error('chebyshev_nodes: d must be a nonnegative integer.');
    end
    

% Allocate output vector
    
    x = zeros(d+1, 1);

% Fill entries using the formula for the roots of T_{d+1}
% i = 0,1,...,d
    
    for i = 0:d
        idx = i + 1;    % Convert i (0-based) to MATLAB index (1-based)

%-------------- Compute x_i = cos( (2*i + 1)*pi / (2*(d+1)) )-----------------------------------------

        numerator   = (2*i + 1) * pi;
        denominator = 2 * (d + 1);
        angle       = numerator / denominator;

        x(idx) = cos(angle);
    end
    
end