function [W, tau] = my_qr_householder(A)
% MY_QR_HOUSEHOLDER  In-place Householder QR with packed storage.
%
% INPUT:
%   A   : n-by-k matrix (with n >= k) with full column rank.
%
% OUTPUT:
%   W   : n-by-k work array that stores:
%           - R in its 'upper-triangle' (rows 1..k, cols 1..k)
%
%           - for each column j, the reflector tail v_j(2:end) in W(j+1:n, j)
%             (v_j is understood to have v_j(1) = 1 because of how we construct v in hh_vec.m)
%
%   tau : k-by-1 vector to store all the scalar multipliers in the Householder transformation 
%               H_j = I - tau(j)*v_j*v_j'
%
% NOTES:
%
% * This routine does NOT form Q. It stores enough data (W, tau) so that
%   later we can apply Q' to any vector b using only basic operations:
%      
% * To extract R afterward, just take its upper triangular block
%       
% * The construction used in hh_vec ensures R has a NONNEGATIVE diagonal.
%
% DEPENDENCIES:
%   - hh_vec.m  (build one reflector with v(1)=1 and tau >= 0)
%
% OUTLINE:
% * At step j, we build a reflector for the active column segment x = A(j:n, j),
%   apply it to the trailing block A(j:n, j:k), then store v_j's tail under
%   the diagonal in column j.
%
% ERROR CHECKS:
% * We require n >= k (tall or square). If n < k, QR with column reflectors
%   does not make sense for least squares as stated here.
%
%=====================================================================================================

    [n, k] = size(A); % Define the size of A

%--------------------Error Message--------------------------------------------------------------------
    if n < k
        error('my_qr_householder: A must have at least as many rows as columns (n >= k).');
    end
%-----------------------------------------------------------------------------------------------------

    W   = A;              % Storage array: will become [R; reflector tails for constructing H_j]
    tau = zeros(k,1);     % store one scalar per reflector/column

    for j = 1:k          % Looping through each column

% -------------- Build a reflector for the active column segment x = W(j:n, j)-----------------------
        
        x = W(j:n, j);           % Define the active column segment of A
        [v, tauj] = hh_vec(x);   % Construct the householder vector used in the Householder 
                                 % Transformation H=I-tauj*(vv') where v(1) = 1 by construction; 
                                 % tauj >= 0
                                 
        tau(j) = tauj;           % store tauj value in the j-th entry of the vect. tau

% ------- Apply Housholder Transformation: H_j = I - tauj * vv' to the trailing block W(j:n, j:k)-----

        if tauj ~= 0                   % If it is 0, the transformation is the identity and we dont 
                                       % have to do anyting
                                       
            X = W(j:n, j:k);                   % Define the working matrix
            w = (v' * X);                      % here v' is acting on each column of X.
            W(j:n, j:k) = X - tauj * (v * w);  % Store the transformation: 
                                               %       HX = [I - tauj*v*v']X = X - tauj(v*w)
                                               % in W.

            
        end

%--------- Store the reflector tail below the diagonal in column j------------------------------------

%      we defined v(1) is implicitly 1; we save only need to store v(2:end)
        
        if j < n
            W(j+1:n, j) = v(2:end);
        end

        % After this step:
        %   - W(j,j:k) stores values of R 
        %   - W(j+1:n,1:j) stores each Householder transformation vector V without the 1 which is 
        %     implicit
        
    end
end