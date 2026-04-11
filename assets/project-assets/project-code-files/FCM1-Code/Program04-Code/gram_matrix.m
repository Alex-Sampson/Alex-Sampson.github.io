function [G, g] = gram_matrix(B, y)
% GRAM_MATRIX  Form the Gram matrix G = B' * B and g = B' * y if y is given.
%
% INPUT:
%   B : n-by-m matrix (e.g., polynomial basis values at n points)
%   y : n-by-1 data vector. If provided, g = B' * y is computed.
%
% OUTPUT:
%   G : m-by-m Gram matrix, G(i,j) = sum_{k=1}^n B(k,i) * B(k,j)
%   g : m-by-1 vector, g(i) = sum_{k=1}^n B(k,i) * y(k)   (only if y is given)
%
% NOTES:
% * For the polynomial regression experiments:
%       - G = B'B is the Gram matrix in the normal equations.
%       - g = B'y is the right-hand side.
% -------------------------------------------------------------------------

    [n, m] = size(B);

    % Initialize G as zeros
    G = zeros(m, m);

  
    for i = 1:m
        for j = 1:m
            s = 0;
            for k = 1:n                    % Compute G = B' * B entry by entry:
                s = s + B(k,i) * B(k,j);   % G(i,j) = sum_{k=1}^n B(k,i) * B(k,j)
            end
            G(i,j) = s;
        end
    end

   
    if nargin < 2
        g = [];                % If y is not provided, we are done
        return;
    end

%---------------------Otherwise, compute g = B' * y-----------------------------
    g = zeros(m, 1);
    for i = 1:m
        s = 0;
        for k = 1:n
            s = s + B(k,i) * y(k);    %g(i) = sum_{k=1}^n B(k,i) * y(k)
        end
        g(i) = s;
    end
end