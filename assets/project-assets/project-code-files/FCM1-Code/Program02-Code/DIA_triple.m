function [DIAG, IOFF] = DIA_triple(n, k, d_main, d_sub, d_super)
% BUILD_DIA_TRIPLE  Create compressed Diagonal storage for a matrix with 3 diagonals:
% offsets -k, 0, +k.  
%
%   DIAG   n-by-3 matrix. Each column is one stored diagonal,
%           padded with zeros so that everything fit's in one
%           n-by-3 array.
%   IOFF   1-by-3 vector of the offsets [ -k, 0, +k ] (same order as columns).
%
% 
%   n        matrix size 
%   k        Offset vallue for the sub/super-diagonals
%   d_main   1-by-n array stores values for the main diagonal  
%   d_sub   1-by-(n-k) array stores values for the sub-diagonal at -k 
%   d_super 1-by-(n-k) array stores values for the super-diagonal at +k 

    % ---- initialize outputs ----
    DIAG = zeros(n, 3);         % Initializes a zero array
    IOFF = [-k, 0, +k];

    % ---- if the diagonals are more than one value, store as column vect's ----
    if length(d_main)   > 1
        d_main   = d_main(:);   
    end 
    if length(d_sub)   > 1
        d_sub   = d_sub(:);   
    end
    if length(d_super) > 1
        d_super = d_super(:); 
    end

  
    if isscalar(d_main)      % if we have just one value on the diag, store it as a repeated n-by-1 array.
        DIAG(:, 2) = d_main * ones(n,1);
    else
        DIAG(:, 2) = d_main; % otherwise, define the diagonal as given.
    end

    %-------------------- fill sub-diagonal (offset -k) -------------------
    % valid rows are i = k+1 through n  (since column = i - k must be >= 1)
    %----------------------------------------------------------------------
    if k > 0
        r1 = (k+1):n;                        % rows that actually have a -k entry
        m  = n - k;                          % number of stored values.
        if isscalar(d_sub)
            DIAG(r1, 1) = d_sub * ones(m,1); % repeat the scalar for all values on sub_diag
        else
            DIAG(r1, 1) = d_sub;             % otherwise define the sub_diag as given.
        end
    else
        % if k == 0, the -k diagonal coincides with main. Do nothing.
    end

    %---------------- fill super-diagonal (offset +k) ---------------------

    % valid rows are i = 1 through n-k  (since column = i + k must be <= n)
    %----------------------------------------------------------------------
    if k > 0
        r2 = 1:(n-k);
        m  = n - k;
        if isscalar(d_super)
            DIAG(r2, 3) = d_super * ones(m,1); % repeat the scalar for all values on the sup_diag
        else
            DIAG(r2, 3) = d_super;             % otherwise define the sup_diag as given.
        end
    else
    % if k == 0, the +k diagonal coincides with main. Do nothing.
    end

end