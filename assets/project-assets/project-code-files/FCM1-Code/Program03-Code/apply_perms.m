function A2 = apply_perms(A, rp, cp)
% APPLY_PERMS  Form A2 = Pr * A * Pc using index vectors.
%   
%     rp : row permutation (length size(A,1)), e.g. [3 1 2]
%     cp : col permutation (length size(A,2)), e.g. [2 3 1]


    if nargin < 2, rp = 1:size(A,1); end
    if nargin < 3, cp = 1:size(A,2); end

    A2 = A(rp, cp);
end
