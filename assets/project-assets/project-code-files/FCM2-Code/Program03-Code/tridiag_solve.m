function x = tridiag_solve(dl, d, du, b)
% tridiag_solve
%   Solve a tridiagonal linear system T*x = b where T is stored as three
%   separate vectors -- exactly the DIA-style storage from FCM1 P2.
%   No full matrix is ever built. O(n) time and O(n) space throughout.
%
%   The system looks like:
%
%   d(1)  du(1)                          x(1)     b(1)
%   dl(1) d(2)  du(2)                    x(2)     b(2)
%         dl(2) d(3)  du(3)         *    x(3)  =  b(3)
%               ...   ...   ...          ...      ...
%                     dl(n-1) d(n)       x(n)     b(n)
%
% Inputs:
%   dl - (n-1) x 1 vector: the sub-diagonal   (dl(i) is in row i+1, col i)
%   d  - (n)   x 1 vector: the main diagonal
%   du - (n-1) x 1 vector: the super-diagonal (du(i) is in row i, col i+1)
%   b  - (n)   x 1 right-hand side vector
%
% Output:
%   x  - (n) x 1 solution vector

    n = length(d);

    % work on copies so we don't mess up the caller's arrays
    d_work = d;
    b_work = b;

    % ---- forward elimination: zero out the sub-diagonal entries ----
    % at each step i we eliminate dl(i-1) from row i using row i-1.
    % this modifies the diagonal d_work(i) and the rhs b_work(i).
    % the super-diagonal du doesn't change during this process because
    % each row only has one entry above the diagonal.
    for i = 2 : n
        factor      = dl(i-1) / d_work(i-1);       % multiplier for this row
        d_work(i)   = d_work(i) - factor * du(i-1); % update diagonal entry
        b_work(i)   = b_work(i) - factor * b_work(i-1); % update rhs
    end

    % ---- back substitution: solve the upper bidiagonal system ----
    % after forward elimination T is upper bidiagonal with d_work on
    % the diagonal and du on the super-diagonal -- solve from the bottom up.
    x = zeros(n, 1);
    x(n) = b_work(n) / d_work(n);   % last unknown is trivial

    for i = n-1 : -1 : 1
        x(i) = (b_work(i) - du(i) * x(i+1)) / d_work(i);
    end

end