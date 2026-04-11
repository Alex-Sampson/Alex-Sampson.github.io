function [A, M, pr, pc, status, rel_factor_err, gamma, rel_sol_err, rel_resid] = driver_set1(n, prob, mode, tol, normtype)
% DRIVER_SET1_ONE  Run ONE Set 1 problem and print/return results.
%   n        : matrix size (default 8)
%   prob     : problem number 1..8 (default 1)
%   mode     : 0=None, 1=Partial (default), 2=Complete
%   tol      : pivot tolerance (default 1e-12)
%   normtype : matrix norm for errors/gamma: 1, Inf, or 'fro' (default 1)
%
% RETURNS:
%   A              original test matrix (kept unchanged for comparison later)
%   M, pr, pc      packed LU and permutations from my_lu(A, mode, tol)
%   status         0=OK, -1=failed pivot (see my_lu)
%   rel_factor_err = ||Pr*A*Pc - (LU)|| / max(1,||A||)
%   gamma          growth factor = || |L||U| || / ||A||
%   rel_sol_err    ||x_hat - x_true||_2 / max(1,||x_true||_2)
%   rel_resid      ||b - A*x_hat||_2 / max(1,||b||_2)

%------------------------------- Defaults --------------------------------
if nargin < 1, n = 8;        end
if nargin < 2, prob = 1;     end
if nargin < 3, mode = 1;     end
if nargin < 4, tol = 1e-12;  end
if nargin < 5, normtype = 1; end

%------------------------------- Build A ---------------------------------
A = zeros(n,n);           % A is the ORIGINAL; do not modify after this point
casename = '(unknown)';

if prob == 1
    casename = 'diag (positive)';
    for i = 1:n, A(i,i) = i; end

elseif prob == 2
    casename = 'antidiag (positive)';
    for i = 1:n, A(i, n-i+1) = i; end

elseif prob == 3
    casename = 'diag + antidiag';
    for i = 1:n
        A(i,i) = i;
        A(i, n-i+1) = n+i;
    end

elseif prob == 4
    casename = 'unit lower, |strict|<1';
    A = eye(n);
    for i = 2:n
        for j = 1:i-1
            A(i,j) = 0.5;
        end
    end

elseif prob == 5
    casename = 'lower (diag>0, strict>1)';
    for i = 1:n
        for j = 1:i
            if i == j, A(i,j) = 1; else, A(i,j) = 2; end
        end
    end

elseif prob == 6
    casename = 'tridiag strictly diag dom';
    for i = 1:n
        A(i,i) = 4;
        if i < n, A(i, i+1) = -1; end
        if i > 1, A(i, i-1) = -1; end
    end

elseif prob == 7
    casename = 'diag=1, below=-1, last col=1';
    for i = 1:n
        A(i,i) = 1;
        A(i,n) = 1;
        for j = 1:i-1, A(i,j) = -1; end
    end

elseif prob == 8
    casename = 'SPD from Ltilde*Ltilde''';
    L = zeros(n,n);
    for i = 1:n
        for j = 1:i
            if i == j, L(i,j) = 1 + 0.5*i/n; else, L(i,j) = 0.2; end
        end
    end
    A = L * L';
else
    error('prob must be an integer in {1,2,...,8}.');
end

%------------------------------ Factorize -------------------------------
[M, pr, pc, status] = my_lu(A, mode, tol);   % A stays unchanged in caller

% Factorization error: ||Pr*A*Pc - (LU)|| / max(1, ||A||)
Aperm = apply_perms(A, pr, pc);
LU    = LUmult_outer(M, n, 0);  % 0 --> LU
denomA = max(1, norm(A, normtype));
rel_factor_err = norm(Aperm - LU, normtype) / denomA;

% Growth factor gamma = || |L||U| || / ||A||
if status < 0
    gamma = NaN;     % factorization aborted
else
    gamma = growth_factor(M, A, normtype);
end

%------------------------------ Solve test ------------------------------
x_true = ones(n,1);
b      = A * x_true;

rel_sol_err = NaN; rel_resid = NaN;
if status >= 0
    x_hat = lu_solve(A, b, mode, tol);
    rel_sol_err = norm(x_hat - x_true, 2) / max(1, norm(x_true, 2));
    rel_resid   = norm(b - A*x_hat,  2) / max(1, norm(b,      2));
end

%------------------------------ Printout --------------------------------
fprintf('-------------------------------------------------------------------------------\n');
fprintf('Set 1  |  n = %d  |  prob = %d (%s)  |  mode = %d  |  tol = %g\n', n, prob, casename, mode, tol);
fprintf('status = %d\n', status);
fprintf('rel factorization error = %.3e\n', rel_factor_err);
fprintf('growth factor gamma     = %.3e\n', gamma);
fprintf('rel solution error      = %.3e\n', rel_sol_err);
fprintf('rel residual            = %.3e\n', rel_resid);

fprintf('\nA (original):\n');    
disp(A);                       
fprintf('\nM (packed LU):\n');  disp(M);
fprintf('pr = ');              disp(pr);
fprintf('pc = ');              disp(pc);
fprintf('-------------------------------------------------------------------------------\n');

end
