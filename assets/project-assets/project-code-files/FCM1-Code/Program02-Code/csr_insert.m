function [AA, JA, NEXT, ok] = csr_insert(i, j, v, AA, JA, IA, NEXT)
% Put (i,j,v) (row, col., val.) at the first free slot of row i if capacity remains.
% Returns ok=true if theres room; ok=false if row i's elbow room is full.
% Note this takes outputs from coo_to_csr_elbow and *updates* the output.

    if NEXT(i) <= IA(i+1)-1       % if we are not at the end of the storage slice
        p = NEXT(i);              % set curser
        AA(p) = v;                % place the value at the curser position
        JA(p) = j;                % store the column index of that position
        NEXT(i) = p + 1;          % advance the curser
        ok = true;                % are we at the end of the storage? No? We can store another value if needed.
    else
        ok = false;               % out of storage room for this row.
    end
end