function y = spmv_modcsr(IA, JA, AA, Diag, x)
% SPMV_MODCSR  Compute y = A*x where A is in Modified CSR:
%   Diag(i) = A(i,i)
%   (IA,JA,AA) = CSR over OFF-diagonal entries only.
% INPUT:
%   IA, JA, AA : off-diagonal CSR (row i uses k = IA(i):IA(i+1)-1)
%   Diag       : diagonal vector (length n)
%   x          : dense vector (length n)
% OUTPUT:
%   y          : dense vector (length n)

     n = length(Diag);       % number of rows
     y = zeros(n,1);         % start result vector as all zeros

    % go through each row
    for i = 1:n
        s = Diag(i) * x(i);           % begin with diagonal entry * x(i)

        k1 = IA(i);                   % first off-diagonal index for row i
        k2 = IA(i+1) - 1;             % last off-diagonal index for row i
        % go through off-diagonal entries of this row
        for k = k1:k2
            j = JA(k);                % column index of this off-diagonal
            v = AA(k);                % value of this off-diagonal
            s = s + v * x(j);         % add its contribution to the sum
        end

        y(i) = s;                     % store the row sum in y(i)
    end
end