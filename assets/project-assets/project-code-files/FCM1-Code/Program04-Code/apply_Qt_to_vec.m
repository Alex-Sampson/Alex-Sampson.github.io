function y = apply_Qt_to_vec(W, tau, b)
    
    
% APPLY_QT_TO_VEC  Compute y = Q' * b using packed Householder array.
%
% INPUT:
%
%   W   : n-by-k work array from my_qr_householder(A).
%        
%        Column j stores:
%           - R(j,j:k) in row j (upper-trapezoid)
%           - the Householder tail v_j(2:end) in rows j+1:n of column j
%
%        The full reflector vector for H_j is v_j = W(j+1:n, j) the j-th col. of W entries j+1:n
%
%   tau : k-by-1 vector of Householder scalars from my_qr_householder.
%
%   b   : n-by-1 right-hand side vector of our linear system of equation.
%
% OUTPUT:

%   y   : The Transformed b vector such that y = Q'*b where
%
%             Q' = H_k H_(k-1) ... H_1
%
%         and each H_j = I - tau(j) * v_j * v_j'.
%
% NOTES:
% * We do NOT build Q explicitly.
% * We start with y = b, then for j = 1..k, apply H_j to the tail y(j:n):
%       y(j:n) <--- y(j:n) - v_j * (tau(j) * (v_j' * y(j:n))).
%
% * This matches the way my_qr_householder applies H_j to the columns of A,
%   so the same sequence produces Q'*b here.
%
% DEPENDENCIES:
%   - my_qr_householder.m   (for how W and tau are produced)
%   - hh_vec.m              (used inside my_qr_householder)

% ====================================================================================================

    [n, k] = size(W);      % Determine the size of W

    y = b(:);              % Make sure that b is a col. vect. store b in y as a work vector
    
%------------------------Error Messages---------------------------------------------------------------

    if length(y) ~= n
        error('apply_Qt_to_vec: length of b must match the number of rows of W.');
    end
    
    if length(tau) ~= k
        error('apply_Qt_to_vec: length of tau must match the number of columns of W.');
    end
    
%-----------------------------------------------------------------------------------------------------

    % Apply each reflector H_j in the same order they were used on A:
    %             y <- H_k ... H_2 H_1 y = Q' * b.
    
    for j = 1:k                    % for each Householder Transformation
        tauj = tau(j);              % define the j-th element of tau
        
        if tauj ~= 0               % if tauj == 0 then the transformation is the identity and we dont       
                                    % need to do anyting    
            v = [1; W(j+1:n, j)];   % Define v from each tail of W's j-th col.
            
            yj = y(j:n);            % Work only on the tail j:n (above that, H_j acts like identity)
            
            s  = v' * yj;           % s = v' * y(j:n) is a scalar

            y(j:n) = yj - v * (tauj * s); % Update y(j:n) <--- y(j:n) - v * (tauj * s)
        end
    end
end