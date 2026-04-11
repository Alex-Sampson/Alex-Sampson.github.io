function [AA, JA, IA, Diag] = coo_to_modcsr(JR, JC, AAc, n)
% COO_TO_MODCSR  Convert in-row-order COO (JR,JC,AAc) to Modified CSR.
% Returns separate arrays:
%   AA   : off diagonal values packed row-by-row
%   JA   : matching column indices for AA
%   IA   : row pointers (len n+1). Row i uses AA/JA indices IA(i) through IA(i+1)-1
%   Diag : diagonal (len n), assumed 0 if missing
%
% Assumptions:
%   - JR, JC, AAc are column vectors of equal length 
%   - Entries are grouped by row (in-row-order)
%   - No duplicate (row,col) pairs
%   - Square matrix size n

    % ensure column vectors
    JR  = JR(:); 
    JC  = JC(:); 
    AAc = AAc(:);

    % number of nonzeros provided
    Nnz = length(AAc);          % TOTAL number of non-zero elements

    % ----- Step 0: initialize diagonal and row off diag arrays -------
    Diag    = zeros(n,1);        % zero array to store diag elements
    cnt_off = zeros(n,1);        % cnt_off(i) stores how many off diagonal elements in row i
    Moff    = 0;                 % TOTAL number of off diagonal entries

    % First pass: split diagonal vs off-diagonal (counts only + fill Diag)
    for k = 1:Nnz                % loops through number of nonzero elements
        i = JR(k);               % row index of non-zero 
        j = JC(k);               % col index of non-zero
        v = AAc(k);              % value of non-zero
        if i == j                % if in the diagonal
            Diag(i) = v;         % store value in Diag (0 stays if missing)
        else %-------------------------- Otherwise --------------------------------------
            cnt_off(i) = cnt_off(i) + 1; % add to the count of the row index
            Moff = Moff + 1;             % add to the count of total off diag elements
        end
    end

    % ---------- Step 1: build IA (row pointers for off diagonals) ----------
    IA    = zeros(n+1,1);               % Initializes our pointer
    IA(1) = 1;                          % Pointer starts as one like before
    for i = 1:n                         % loops through the rest of the pointer indices
        IA(i+1) = IA(i) + cnt_off(i);   % assigns the next slot to be the value of the previous slot + num of 
    end                                 % non-zero off diag entries.
    
    % Now row i occupies AA/JA indices IA(i) through  IA(i+1)-1.
    % Note: empty row i if IA(i) = IA(i+1).

    % ---------- Step 2: initialize off diagonal storage and row cursors -----
    AA   = zeros(Moff,1);               % zero aray of length Moff for off diag val's
    JA   = zeros(Moff,1);               % zero aray of length Moff for col index
    next = IA(1:n);                     % next(i) = next free slot in row i

    % ---------- Step 3: second pass to write off diagonals ----------------
    for k = 1:Nnz              % SAME LOOPING AS STEP 0 - First pass.
        i = JR(k);   
        j = JC(k);
        v = AAc(k);
        if i ~= j              % if i =/= j, i.e. off diagonal
            p    = next(i);    % position to write this off diagonal
            AA(p) = v;         % value
            JA(p) = j;         % column index
            next(i) = p + 1;   % advance row i cursor
        end
    end
end