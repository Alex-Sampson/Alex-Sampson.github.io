function [x, R, c, d, W, tau] = lsq_householder_incremental(A, b, tol)
% LSQ_HOUSEHOLDER_INCREMENTAL
%
% Solve the least-squares problem
%
%        minimize  || b - A x ||_2
%
% using Householder QR, but apply each reflector H_j to b individually
% during the factorization, instead of calling apply_Qt_to_vec afterwards.
%
% INPUT:
%   A   : n-by-k matrix (tall or square, with n >= k, full column rank)
%   b   : n-by-1 right-hand side vector
%   tol : tolerance for back_sub 
%
% OUTPUT:
%   x   : least-squares minimizer
%   R   : k-by-k upper-triangular factor (top block of W)
%   c   : first k entries of y = Q' * b (computed incrementally)
%   d   : last  n-k entries of y = Q' * b
%   W   : n-by-k work array containing R and packed Householder vectors
%   tau : k-by-1 vector of Householder scalars
%
% NOTES:
% * This routine uses the same reflector construction as my_qr_householder,
%   but we explicitly update b at each step:
%        y <--- H_j y
%   so that by the end, y = Q' b.
%
% * This is mainly for comparison with LSQ_HOUSEHOLDER to show that
%   factorization and "incremental application of transformations produce the same c, d, and x.
% ====================================================================================================

    if nargin < 3
        tol = 1e-12;  % set default tol
    end

    [n, k] = size(A); % define the size of A
    
%-------------------------------Error Message---------------------------------------------------------
    
    if n < k
        error('lsq_householder_incremental: A must have n >= k.');
    end

%-----------------------------------------------------------------------------------------------------

%----------------------Initialize work array and tau, and copy b into y-------------------------------
    W = A;
    tau = zeros(k,1);
    y = b(:);   % y will be updated to H_k ... H_1 b = Q' b
%-----------------------------------------------------------------------------------------------------

    for j = 1:k

% --------------------Build reflector for the active column segment W(j:n, j)-------------------------
        x_col = W(j:n, j);
        [v, tauj] = hh_vec(x_col);   % v(1)=1, tauj >= 0
        tau(j) = tauj;

%------------------------Apply H_j to the trailing block W(j:n, j:k)----------------------------------
        if tauj ~= 0
            X = W(j:n, j:k);         % active block of A
            w_row = v' * X;          % 1-by-(k-j+1)
            W(j:n, j:k) = X - tauj * (v * w_row);
        end

%-------------------- Apply the same H_j to y (incremental update of Q' b)----------------------------
        if tauj ~= 0
            y_tail = y(j:n);         % tail y(j:n)
            s = v' * y_tail;         % scalar
            y(j:n) = y_tail - tauj * (v * s);
        end

%-------------------Store reflector tail below the diagonal in column j-------------------------------
        if j < n
            W(j+1:n, j) = v(2:end);
        end
    end

% -------------------- Extract R (top k-by-k block) --------------------------------------------------
    R=zeros(k,k);
    
    for i=1:k
        for j=1:k
           R(i,j) = W(i, j);
        end
    end

% -------------------- Split y = [c; d] --------------------------------------------------------------
    c = y(1:k);
    d = y(k+1:n);

% -------------------- Solve R x = c by back-substitution --------------------------------------------
    x = back_sub(R, c, tol);
end