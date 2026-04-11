function [VAL, COL, K] = coo_to_ellpack(JR, JC, AAc, n)
% COO_TO_ELLPACK  Convert in-row-order COO (JR,JC,AAc) to ELLPACK/ITPACK.
%
% Returns:
%   VAL : n-by-K array of values (each row padded with 0's to length K)
%   COL : n-by-K array of column indices (padded with 0's)
%   K   : max number of stored entries in any row (row capacity)
%
% Assumptions:
%   - JR, JC, AAc are column vectors of equal length 
%   - Entries are grouped by row (in-row-order)
%   - No duplicate (row,col) pairs
%   - n = number of rows (matrix is square, but ELLPACK does not require it)

    % --- Step 0: ensure column vector format ------------------------
    JR  = JR(:);   % make sure JR is a column vector
    JC  = JC(:);   % make sure JC is a column vector
    AAc = AAc(:);  % make sure AAc is a column vector

    % --- Step 1: initialize counters -------------------------------
    Nnz = length(AAc);   % total number of nonzero entries in COO
    cnt = zeros(n,1);    % cnt(i) will hold number of entries in row i
    K   = 0;             % K = maximum row length (row capacity)

    % --- Step 2: first pass through COO to count per row and find K -
    for k = 1:Nnz             % loop over every nonzero in COO
        i = JR(k);            % row index of current entry
        cnt(i) = cnt(i) + 1;  % bump the count for that row
        if cnt(i) > K         % check if this row is now the largest so far
            K = cnt(i);       % update K if needed
        end
    end
    % After this loop:
    %   cnt(i) = how many entries row i has
    %   K      = max row length (capacity of each row in ELLPACK)

    % --- Step 3: Initialize storage locations for VAL and COL ------------------
    VAL  = zeros(n, K);    % each row padded with 0 up to length K
    COL  = zeros(n, K);    % each row's column indices, padded with 0
    next = ones(n,1);      % next(i) = next free slot in row i (starts at 1)

    % --- Step 4: second pass to actually store VAL and COL ----
    for k = 1:Nnz
        i = JR(k);         % row index of this entry
        j = JC(k);         % column index of this entry
        v = AAc(k);        % value of this entry

        p = next(i);       % which slot in row i we should fill
        VAL(i,p) = v;      % store the value into row i, slot p
        COL(i,p) = j;      % store the column index into row i, slot p
        next(i)   = p + 1; % advance the cursor for row i to next free slot
    end
    % After this loop:
    %   Each row i of VAL contains its nonzeros (in-row-order), padded with 0's
    %   Each row i of COL contains the column indices, padded with 0's
end
