function x = back_sub(M, y, tol)
% Solve U*x = y where U is upper triangular.

if nargin < 3, tol = 1e-12; end

[n,m] = size(M);

if n ~= m, error('M must be square.'); end

y = y(:); % insure y is a column vector

if length(y) ~= n, error('y has wrong length.'); end

x = zeros(n,1); % initialize x

for i = n:-1:1                                          % loops BACKWARDS from n to 1 with integer steps since we will work bottom up on U
    s = 0;                                              % firs step requires no subtraction
    for j = i+1:n                                       % for row i, loop through each off diag col element
        s = s + M(i,j)*x(j);                            % update s to account for all col entries in row i
    end
    
    piv = M(i,i);                                       % or divisor for getting final solution. Pivot (Diag) entry of U.
    
    if abs(piv) < tol
        error('Zero/near-zero pivot at row %d.', i); 
    end
    
    x(i) = (y(i) - s)/piv;                              % solve for x(i)
    
end
end