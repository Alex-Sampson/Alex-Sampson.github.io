function y = Sym_CSR_spmv(IA, JA, AA, x)
% y = B*x for a symmetric matrix B stored in CSR
% that computes using ONLY the lower triangle & diagonal (j <= i).
%
%   - IA: row pointer, IA(i): start index of row i in JA/AA 
%   - JA: column indices for stored entries (only j <= i)
%   - AA: values for those entries
%   - Row i occupies positions p = IA(i) through IA(i+1)-1 in JA/AA


    n = length(IA) - 1;  % Our row pointer is always one element longer than the size of B.
    y = zeros(n,1);     % initialize y as a zero array.

    for i = 1:n
        s_i = 0;                        % sum to compute y(i)
        for p = IA(i):IA(i+1)-1         % for our slice of row i
            j = JA(p);                  % store column index
            a = AA(p);                  % store value
%==============================================================================
% Strictly lower diagonal entries will contribute two entries, due to symmetry
%==============================================================================
            if j < i                       % below the diagonal
                s_i   = s_i   + a * x(j);  % to compute y(i)
                y(j) = y(j) + a * x(i);    % symmetric update to y(j)
            else                           % otherwise, if a diagonal entry
                s_i   = s_i   + a * x(i);  % sum is computed as usual, no symmetric update
            end
        end
        y(i) = y(i) + s_i;                 % compute y(i)
    end
end