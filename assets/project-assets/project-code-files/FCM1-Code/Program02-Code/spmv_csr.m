function y = spmv_csr(IA, JA, AA, x)
% SPMV_CSR  Compute y = A*x where A is in standard CSR format.
% INPUT:
%   IA : row pointer (length n+1). Row i uses k = IA(i):IA(i+1)-1
%   JA : column indices for AA
%   AA : nonzero values (aligned with JA)
%   x  : dense vector (length n)
% OUTPUT:
%   y  : dense vector (length n)

    n = length(IA) - 1;      % number of rows in the matrix
    y = zeros(n,1);          % start result vector as all zeros

    % go through each row of the matrix
    for i = 1:n
        s = 0;                        % this will hold the sum for row i
        k1 = IA(i);                   % first position in AA/JA for row i
        k2 = IA(i+1) - 1;             % last position in AA/JA for row i
        % go through each stored entry in row i
        for k = k1:k2
            j = JA(k);                % column index of this entry
            v = AA(k);                % value of this entry
            s = s + v * x(j);         % add its contribution to the row sum
        end
        y(i) = s;                     % store the final sum in y(i)
    end
end