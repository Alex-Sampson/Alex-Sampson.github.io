function [alpha, breaks] = spline2_coeffs(x_nodes, f_nodes, bc_type, bc_vals)
% spline2_coeffs
%   Build an interpolatory cubic spline using the B-spline basis on a
%   UNIFORM mesh. Solves the (n+3) x (n+3) linear system A*alpha = b
%   for the B-spline coefficients alpha_{-1}, alpha_0, ..., alpha_{n+1}.
%
%   The spline is written as:
%       s(x) = sum_{i=-1}^{n+1} alpha_i * B_i(x)
%
%   where B_i is the cubic B-spline centered at x_i (Definition 3.1,
%   Set 14 notes). We restrict to a uniform mesh as permitted by the
%   project outline.
%
%   Boundary conditions:
%   The system has n+1 interpolation equations (one per knot) plus 2 BC
%   equations. Both BC types below are supported:
%
%   'clamped' : s'(x_0) = bc_vals(1),  s'(x_n) = bc_vals(2)
%               -- the first derivative rows use B_i'(x_j) from Lemma 3.2
%
%   'natural' : s''(x_0) = 0,  s''(x_n) = 0
%               -- second derivative rows from differentiating B_i twice
%               -- B_{-1}''(x_0)=6/h^2, B_0''(x_0)=-12/h^2, B_1''(x_0)=6/h^2
%
% Inputs:
%   x_nodes  - (n+1) x 1 vector of UNIFORMLY spaced knot locations
%   f_nodes  - (n+1) x 1 vector of function values f(x_i)
%   bc_type  - 'clamped' or 'natural'
%   bc_vals  - [left_val, right_val] for clamped; ignored for natural
%
% Outputs:
%   alpha  - (n+3) x 1 vector of B-spline coefficients
%   breaks - same as x_nodes, returned for spline2_eval

    n = length(x_nodes) - 1;      % number of subintervals
    h = x_nodes(2) - x_nodes(1);  % uniform spacing

    N = n + 3;   % total unknowns: alpha_{-1} .. alpha_{n+1}

    A = zeros(N, N);
    b = zeros(N, 1);

    % ---- index convention ----
    % alpha_{-1} -> position 1
    % alpha_0    -> position 2
    % alpha_i    -> position i+2
    % alpha_{n+1}-> position N

    % ---- first row: left BC ----
    if strcmp(bc_type, 'clamped')
        % s'(x_0) = bc_vals(1)
        % B_{-1}'(x_0) = -3/h,  B_0'(x_0) = 0,  B_1'(x_0) = 3/h
        A(1, 1) = -3/h;
        A(1, 3) =  3/h;
        b(1)    = bc_vals(1);

    elseif strcmp(bc_type, 'natural')
        % s''(x_0) = 0
        % B_{-1}''(x_0) = 6/h^2,  B_0''(x_0) = -12/h^2,  B_1''(x_0) = 6/h^2
        A(1, 1) =  6 / h^2;
        A(1, 2) = -12 / h^2;
        A(1, 3) =  6 / h^2;
        b(1)    = 0;

    else
        error('bc_type must be ''clamped'' or ''natural''');
    end

    % ---- rows 2 through N-1: interpolation conditions s(x_j) = f_j ----
    % at knot x_j: B_{j-1}(x_j)=1, B_j(x_j)=4, B_{j+1}(x_j)=1
    for j = 0 : n
        row        = j + 2;
        A(row, j+1) = 1;
        A(row, j+2) = 4;
        A(row, j+3) = 1;
        b(row)      = f_nodes(j+1);
    end

    % ---- last row: right BC ----
    if strcmp(bc_type, 'clamped')
        % s'(x_n) = bc_vals(2)
        % B_{n-1}'(x_n) = -3/h,  B_n'(x_n) = 0,  B_{n+1}'(x_n) = 3/h
        A(N, N-2) = -3/h;
        A(N, N)   =  3/h;
        b(N)      = bc_vals(2);

    elseif strcmp(bc_type, 'natural')
        % s''(x_n) = 0
        % B_{n-1}''(x_n) = 6/h^2, B_n''(x_n) = -12/h^2, B_{n+1}''(x_n) = 6/h^2
        A(N, N-2) =  6 / h^2;
        A(N, N-1) = -12 / h^2;
        A(N, N)   =  6 / h^2;
        b(N)      = 0;
    end

    % ---- forward elimination (no pivoting -- A is diagonally dominant) ----
    for k = 1 : N-1
        for i = k+1 : N
            if A(i, k) ~= 0
                factor   = A(i, k) / A(k, k);
                for j = k : N
                    A(i, j) = A(i, j) - factor * A(k, j);
                end
                b(i) = b(i) - factor * b(k);
            end
        end
    end

    % ---- back substitution ----
    alpha = zeros(N, 1);
    for i = N : -1 : 1
        alpha(i) = b(i);
        for j = i+1 : N
            alpha(i) = alpha(i) - A(i, j) * alpha(j);
        end
        alpha(i) = alpha(i) / A(i, i);
    end

    breaks = x_nodes;

end