function y = forward_unit(M, b)
% Solve L*y = b where L is unit lower triangular.

[n,m] = size(M);
if n ~= m, error('M must be square.'); end

b = b(:); % insure b is a column vector

if length(b) ~= n, error('b has wrong length.'); end


y = zeros(n,1);               % initialize y as a zero vect.
for i = 1:n                   % loop through all rows
    s = 0;                    % firs step requires no subtraction
    for j = 1:i-1             % loop col elements befor diag
        s = s + M(i,j)*y(j);  % update s to account for all col entries in row i
    end
    y(i) = b(i) - s;          % solve for i-th element of i
end
end

