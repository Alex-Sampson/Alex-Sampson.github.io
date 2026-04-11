function gamma = growth_factor(M, A0, normtype)

% This function calculates the growth factor of a matrix A after it has been factored into its LU decomposition.
% M is the packed LU matrix, and A0 is the original matrix A that was factored.

% Check if the dimensions of M and A0 are compatible for multiplication
if size(M, 1) ~= size(A0, 1)
    error('The number of rows in M must match the number of rows in A0.'); % Raise an error if dimensions do not match
end

    if nargin < 3
        normtype = 1;
    end

    n = size(M,1);
    AbsLU = LUmult_outer(M, n, 1);   % 1 → |L||U|
    gamma = norm(AbsLU, normtype) / norm(A0, normtype);

   
end
