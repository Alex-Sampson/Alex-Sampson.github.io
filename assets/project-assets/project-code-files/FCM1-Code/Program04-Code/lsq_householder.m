function [x, R, c, d, W, tau] = lsq_householder(A, b, tol)
% LSQ_HOUSEHOLDER  Solve the least-squares problem:
%
%        minimize  || b - A x ||_2
%
% using Householder QR:

%        A = Q R,   R upper-triangular with NONNEGATIVE diagonal.
%        Q' b = [c; d]
%
%     ||b-Ax||=||Q'b-Q'Ax||=||[c;d]-[R_1x;0]|| 
%
%-----Least Squares Problem--------------------------> minimize:||c-R_1x||
%
%    R_1 is non-singular so the unique solution is: x = R_1^{-1} c
%       
%
% INPUT:
%   A   : n-by-k matrix (tall or square, with n >= k)
%   b   : n-by-1 vector
%   tol : small tolerance for back_sub (optional)
%
% OUTPUT:
%   x   : least-squares minimizer
%   R   : k-by-k upper-triangular R factor (extracted from W)
%   c   : first k entries of y = Q' b
%   d   : last  n-k entries of y = Q' b   (the residual component)
%   W,tau : packed QR data from my_qr_householder
%
% NOTES:
% * This routine does NOT form Q explicitly; it uses apply_Qt_to_vec.
% * The diagonal of R is guaranteed nonnegative because hh_vec chooses
%   the sign of the reflector.
%
% ====================================================================================================

    if nargin < 3
        tol = 1e-12; % Set default tollerance
    end

    % Step 1: Factor A using our in-place Householder QR
    
    [W, tau] = my_qr_householder(A);

    % Step 2: Form y = Q' * b using packed reflectors
    
    y = apply_Qt_to_vec(W, tau, b);

    % Step 3: Split y = [c; d], where c has length k
    
    [n, k] = size(A); % Determine the size of A
    
    c = y(1:k);       % Extract c from y=Q'b
    d = y(k+1:n);     % Extract d from y=Q'b

    % Step 4: Extract R from the top k-by-k block of W
    R=zeros(k,k);
    
    for i=1:k
        for j=1:k
           R(i,j) = W(i, j);
        end
    end

    % Step 5: Solve R x = c using our upper-triangular solver from poject 03.
    
    x = back_sub(R, c, tol);
end