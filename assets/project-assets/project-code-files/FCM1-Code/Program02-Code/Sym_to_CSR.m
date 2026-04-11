function [IA, JA, AA] = Sym_to_CSR(B)
% Converts a Symmetric n-by-n matrix to CSR storage format.
% Store ONLY the diagonal and lower triangle (j <= i) of a symmetric B

    n  = size(B,1);      % Define problem size
    IA = zeros(n+1,1);   % Initialize row pointer

    % ---- non-zero count per row (lower & diag) ----
    for i = 1:n
        cnt = 0;                 % start count at zero
        for j = 1:i              % loop through lower half entries
            if B(i,j) ~= 0       % if entry is non-zero
                cnt = cnt + 1;   % increas count by 1
            end
        end
        IA(i) = cnt;             % store non-zero count for lower half of B
    end

    % ---- Construct row pointer ----
    position = 1;                   % start of JA/AA
    for i = 1:n
        cnt      = IA(i);           % old count for row i
        IA(i)    = position;        % start index for row i
        position = position + cnt;  % define starting position for next row
    end
    IA(n+1) = position;             % one past the end

    % total kept nonzeros in the lower half
    nnzL = IA(n+1) - 1;            

    % ---- fill JA/AA using row curser ----
    JA   = zeros(nnzL,1);        % Initialize column index array as zeros
    AA   = zeros(nnzL,1);        % Initialize value array as zeros
    next = IA;                   % next(i) = next free slot for row i

    for i = 1:n
        for j = 1:i                  % looping through lower half of B
            aij = B(i,j);            % store value at row i, col. j
            if aij ~= 0              % if that value is non-zero
                p    = next(i);      % define curser
                JA(p) = j;           % store column index of element
                AA(p) = aij;         % store vlaue of element
                next(i) = p + 1;     % advance the curser one step.
            end
        end
    end
end