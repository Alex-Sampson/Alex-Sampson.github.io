function x = lu_solve(A, b, mode, tol)
    
% Solve A x = b with packed LU and index arrays.
%        mode: 
%            0 --> None 
%            1 --> Partial --> Default
%            2 --> Complete

%------------------------ERROR'S/DEFAULTS---------------------------------

if nargin < 3, mode = 1; end 

if nargin < 4, tol  = 1e-12; end

[n,m] = size(A);

if n ~= m
    error('A must be square.');
end

b = b(:);

if length(b) ~= n
    error('b must have length n.');
end

%------------------------------------------------------------------------
[M,pr,pc] = my_lu(A, mode, tol);

bhat = b(pr);                          % apply row perm (P*b)

y = forward_unit(M, bhat);             % forward on L (packed in M)

z = back_sub(M, y, tol);               % back on U (packed in M)

if mode == 2                           % undo column perm: x = (Pc)*z
    x = zeros(n,1);
    for i = 1:n
        x(pc(i)) = z(i);
    end
    
else
    x = z;
end
end