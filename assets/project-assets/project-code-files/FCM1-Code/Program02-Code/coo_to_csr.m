function [AA, JA, IA] = coo_to_csr(JR, JC, AAc, n) 
%==========================================================================
% I wrote this function not thinking we were going to be given the AAc, JR
% and JC data in-row-order. It should be robust to out of order data. So we
% will use it in TASK2 as well.
%==========================================================================

% COO_TO_CSR  Convert in-row-order COO (JR,JC,AAc) to CSR (IA,JA,AA).
%  (Task 1):
%   - JR, JC, AAc are column vectors of equal length.
%   - Entries are in-row-order (rows grouped).
%   - No duplicate (row,col) pairs.
%   - 1-based indices; n = number of rows.
%
% Returns:
%   AA(k): k-th stored nonzero value in row order
%   JA(k): column index for AA(k)
%   IA(i): starting pointer for row i (IA has length n+1)

    % ensure column vectors
    JR  = JR(:); 
    JC  = JC(:); 
    AAc = AAc(:);

    % number of nonzeros in our Matrix
    Nnz = length(AAc);

    %----- Step 1: count nnz per row
    cnt = zeros(n,1);
    for k = 1:Nnz                     % loop over all nonzeros
        cnt(JR(k)) = cnt(JR(k)) + 1;  % bump the count for that row
    end

    %----- Step 2: build IA (row pointer), 
    IA    = zeros(n+1,1);
    IA(1) = 1;                         % row 1 starts at AA/JA index 1
    for i = 1:n
        IA(i+1) = IA(i) + cnt(i);      % cumulative sum of row lengths
    end
    
    % Now row i occupies AA/JA indices IA(i) through IA(i+1)-1.

    %----- Step 3: Initialize AA/JA and row cursor
    AA   = zeros(Nnz,1);               
    JA   = zeros(Nnz,1);
    next = IA(1:n);                    % next(i) = next free slot for row i

    %----- Step 4: Build AA, JA and move next.
    for k = 1:Nnz
        i = JR(k);       % row index of this nonzero
        j = JC(k);       % column index of this nonzero
        v = AAc(k);      % value of this nonzero

        p = next(i);     % next free position inside row i's slice
        AA(p) = v;       % store the value
        JA(p) = j;       % store the COLUMN index 
        next(i) = p + 1; % advance the row i cursor to the next slot
    end
end