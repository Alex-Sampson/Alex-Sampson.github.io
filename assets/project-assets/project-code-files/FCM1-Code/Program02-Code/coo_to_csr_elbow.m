function [AA, JA, IA, NEXT] = coo_to_csr_elbow(JR, JC, AAc, n, K)
% COO_TO_CSR_ELBOW  Unordered COO -> CSR with per-row elbow room of size K.
%   JR, JC, AAc : column vectors 
%   n           : size of n-by-n matrix
%   K           : constant padding size for all rows
%
% Returns:
%   IA(i)    : start of row i's capacity slice (length n+1)
%   AA, JA   : sized to total capacity
%   NEXT(i)  : first free slot in row i (so used = IA(i):NEXT(i)-1)

    if nargin < 5, K = 4; end % If I dont include a specific padding value we default to 4.

   
    Nnz = numel(AAc); % counts nonzero elements

    % ---- count nnz per row ----
    cnt = zeros(n,1);                     % Initialized as zero array
    for k = 1:Nnz                         % loop through total number of nonzero elements
        cnt(JR(k)) = cnt(JR(k)) + 1;      % every time a row index comes up, add a count to that slot in cnt.
    end

    % ---- capacities and IA (capacity pointers) ----
    cap = cnt + K;                        % the capacity of each row is the total nnz + whatever elbow room we add
    IA  = zeros(n+1,1);                   % Remember: IA(i) start of row i, IA(i+1) one past the end or row i.
    IA(1) = 1;                            % IA is the row pointer but tells us nonzeros + storage for each row
    for i = 1:n
        IA(i+1) = IA(i) + cap(i);         % Defines the start of the next row (after accounting for storage capacity)
    end
    C = IA(end) - 1;                      % total capacity of the matrix should be equal to nnz+n*K !IF! K is const.
                                          % across all rows

    % ---- storage ----
    AA   = zeros(C,1);          % Initialize Storage location for matrix values
    JA   = zeros(C,1);          % Column index of each entry
    NEXT = IA(1:n);             % Initialized to IA but will be updated to point to first available storage slot.

    % ---- fill original nnz at the *front* of each row slice ----
    for k = 1:Nnz               % Loop through all nonzero entries
        i = JR(k);              % row index of the k-th nonzero entry
        p = NEXT(i);            % the first free slot in row i's slice (a curser)
        AA(p) = AAc(k);         % stores the k-th nonzero entry in that free slot
        JA(p) = JC(k);          % Assigns the column index of that entry
        NEXT(i) = p + 1;        % moves the curser forward one step.
    end
end
