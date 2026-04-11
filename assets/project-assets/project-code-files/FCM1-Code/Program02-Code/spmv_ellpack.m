function y = spmv_ellpack(VAL, COL, K, x)
% SPMV_ELLPACK  Compute y = A*x where A is in ELLPACK/ITPACK format.
% INPUT:
%   VAL : n-by-K values (row-padded with 0)
%   COL : n-by-K column indices (0 in padded slots)
%   K   : row capacity (max nonzeros per row)
%   x   : dense vector (length n)
% OUTPUT:
%   y   : dense vector (length n)

    n = size(VAL, 1);       % number of rows
    y = zeros(n,1);         % start result vector as all zeros

    % go through each row
    for i = 1:n
        s = 0;                        % this will hold the sum for row i
        % go through all K slots in this row
        for p = 1:K
            j = COL(i,p);             % column index at this slot
            if j == 0                 % if column index is 0, it is just padding
                % do nothing, skip this slot
            else
                v = VAL(i,p);         % value stored in this slot
                s = s + v * x(j);     % add its contribution to the row sum
            end
        end
        y(i) = s;                     % store the row sum in y(i)
    end
end