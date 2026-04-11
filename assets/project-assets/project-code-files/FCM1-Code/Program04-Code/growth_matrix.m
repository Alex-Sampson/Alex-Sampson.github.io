function A = growth_matrix(n)
    
% GROWTH_MATRIX  Construct the "growth" test matrix used in Program 3.
%
% For an n-by-n matrix A:
%   - A(i,i)   = 1       (diagonal entries)
%   - A(i,j)   = -1      for j < i (strictly below diagonal)
%   - A(i,n)   = 1       (last column all ones)
%
% This is the matrix that produced large growth factors for LU in Program 3.
% In Program 4, we use it to compare the size of entries in the R factor
% from Householder QR with the U factor from LU.

    A = zeros(n,n);

    for i = 1:n
        A(i,i) = 1;          % diagonal
        A(i,n) = 1;          % last column
        for j = 1:i-1
            A(i,j) = -1;     % below diagonal
        end
    end
end