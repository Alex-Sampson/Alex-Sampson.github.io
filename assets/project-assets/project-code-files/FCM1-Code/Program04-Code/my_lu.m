function [M,pr,pc,status] = my_lu(M, mode, tol)

% Gaussian elimination with complete pivoting
% Mode selection for choice of pivoting:
%    0 --> None
%    1 --> Partial
%    2 --> Complete


%--------------------------Error Message/Defaults ---------------------------------------
if nargin < 3
    tol = 1e-12; % if tolerance not specified, use default
end

if nargin < 2
    mode=1;
end
   
status=0; % status = 0 --> Ok, = 1 --> all pivot candidates < tol, = -1 --> fail (pivot < tol)

[n,m] = size(M);
if n ~= m
    error('M must be square.');
end

%-------------------------------------------------------------------------------

pr=1:n; pc=1:n; %initialize row/col pivot integer arrays


for k = 1:n-1 % looping through pivot rows (last row not needed)
    
    if mode==0

%========================NO PIVOT==============================================================

                                     
    piv = M(k,k);    % define a pivot value

    if abs(piv) < tol                                
        warning('Pivot below tol at step %d (no pivoting).', k);   
        status = 1;                                    % set status value
    end


%========================END NO PIVOT==============================================================









    elseif mode==1
%========================PARTIAL PIVOT==============================================================

% ------------------------choose pivot row-----------------------------
    p = k;                           % set pivot position
    maxval = abs(M(k,k));            % set max value to the diagonal element, our potential pivot
    for r = k+1:n                    % loop through rows below pivot entry
        if abs(M(r,k)) > maxval      % if any of those elements are greater than max val
            maxval = abs(M(r,k));    % update max val
            p = r;                   % record row index for that max val
        end
    end

    if maxval < tol                                          % all candidates in column k are < tol
        warning('All candidate pivots in column %d below tol at step %d (partial).', k, k);  
        status = 1;                                          % set status value
    end

%------------ swap rows k and p in M and P (also moves prior L multipliers) --------------

    if p ~= k                                         % if the max entry is not in the pivot spot
        tmp = M(k,:); M(k,:) = M(p,:); M(p,:) = tmp;  % store row k of M temporarily, update row k to be row p, update M
        tmp = pr(k); pr(k) = pr(p); pr(p) = tmp;  % store element k of pr temporarily, update element k to element p of pr, update pr array.
    end

    piv = M(k,k);                                     % now that we know the pivot is the maximal element, store it.
    

%========================END PARTIAL PIVOT===========================================================







    elseif mode==2
%========================COMPLETE PIVOT==============================================================

% ------------------------- Search for max entry ---------------------------------------

    p = k; q = k; maxval = abs(M(k,k));     % for each row initialize p,q and maxval
    for i = k:n                             % loop through rows below pivot entry
        for j = k:n                         % loop through column entries of each of the row entries below pivot (i.e. loop through submatrix)
            v = abs(M(i,j));                % store max value found in this submatrix
            if v > maxval                   % if this max val is greater than initialized value
                maxval = v; p = i; q = j;   % update the max val and store its row/col index
            end
        end
    end

    if maxval < tol                                           % all candidates in trailing submatrix are < tol
        warning('All candidate pivots in trailing submatrix below tol at step %d (complete).', k);  
        status = 1;                                           % set status value
    end

    % ----------------------------- swap rows k and p ------------------------------------------
    if p ~= k                                           % if it wasnt in the pivot position
        tmp = M(k,:); M(k,:) = M(p,:); M(p,:) = tmp;    % store pivot row, swap pivot row with row with max element, update row p to row k's vals.
        tmp = pr(k); pr(k) = pr(p); pr(p) = tmp;    % temp store row k of pivot matrix P , update row k to row p's vals, update row p to row k's vals
    end

    % ----------------------------- swap cols k and q ------------------------------------------
    if q ~= k                                           % if it wasnt in the pivot position
        tmp = M(:,k); M(:,k) = M(:,q); M(:,q) = tmp;    % store the pivot col, update col k to col q's vals, update col q to col k's vals.
        tmp = pc(k); pc(k) = pc(q); pc(q) = tmp;    % temp store element k of col pivot integer array pc, update element k to element q's val, update element q to element k's val.
    end

    piv = M(k,k);                                       % define the pivot value
    

%========================END COMPLETE PIVOT==========================================================
    else
        error('Must choose a pivot mode: 0=None, 1=Partial, 2=Complete');
    end
    
    % ---------------------------------- Error Message ---------------------------------------
    
    if abs(piv) < tol
        status = -1;                % Fail: (no acceptable pivot)
        return                      
    end

    % ----------------------------- multipliers below pivot ------------------------------------------
    for i = k+1:n               % for each col k
        M(i,k) = M(i,k) / piv;  % define the multipliers below that col.
    end

    % ----------------------------- Update M by doing the elimination -----------------------------
    for i = k+1:n
        M(i,k+1:n) = M(i,k+1:n) - M(i,k) * M(k,k+1:n); % subtract the multiple of the pivot row from each row below pivot row 
    end
end
end