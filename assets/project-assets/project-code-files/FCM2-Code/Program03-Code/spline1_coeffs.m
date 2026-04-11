function [s2, breaks, f_vals] = spline1_coeffs(x_nodes, f_nodes, bc_type, bc_vals)
% spline1_coeffs
%   Build an interpolatory cubic spline using the s_i'' (second-derivative)
%   parameterization. Assembles the tridiagonal system T*s2 = d for the
%   unknown second derivatives at each knot and solves it using tridiag_solve.
%
%   T is stored as three vectors (sub-diagonal, main diagonal, super-diagonal)
%   in the DIA style from FCM1 P2 -- no full matrix is ever built.
%   Supports nonuniform meshes throughout.
%
% Inputs:
%   x_nodes  - (n+1) x 1 vector of strictly increasing knot locations
%   f_nodes  - (n+1) x 1 vector of function values f(x_i)
%   bc_type  - string specifying the boundary condition type:
%                'natural'  : s''(a) = 0 and s''(b) = 0
%                'clamped'  : s'(a) and s'(b) are prescribed (see bc_vals)
%   bc_vals  - 1 x 2 vector [left_val, right_val].
%                For 'natural'  : ignored, pass [] if you like
%                For 'clamped'  : bc_vals = [s'(a), s'(b)]
%
% Outputs:
%   s2       - (n+1) x 1 vector of second derivatives s''(x_i) at each knot
%   breaks   - same as x_nodes, returned so spline1_eval doesn't need to
%              carry x_nodes separately
%   f_vals   - same as f_nodes, returned for the same reason

    n = length(x_nodes) - 1;   % number of subintervals

    % ---- subinterval widths h_i = x_i - x_{i-1}, for i = 1..n ----
    h = zeros(n, 1);
    for i = 1 : n
        h(i) = x_nodes(i+1) - x_nodes(i);
    end

    % ---- right-hand side d_i for interior knots i = 1..n-1 ----
    % d_i = 6 * f[x_{i-1}, x_i, x_{i+1}]  (scaled second divided difference)
    % = 6/(h_i + h_{i+1}) * ((f_{i+1}-f_i)/h_{i+1} - (f_i-f_{i-1})/h_i)
    d_int = zeros(n-1, 1);
    for i = 1 : n-1
        left_diff  = (f_nodes(i+1) - f_nodes(i))   / h(i);
        right_diff = (f_nodes(i+2) - f_nodes(i+1)) / h(i+1);
        d_int(i)   = 6 * (right_diff - left_diff) / (h(i) + h(i+1));
    end

    % ---- off-diagonal weights for the interior rows ----
    % mu_i  = h_i     / (h_i + h_{i+1})   multiplies s2_{i-1}  (sub-diagonal)
    % lam_i = h_{i+1} / (h_i + h_{i+1})   multiplies s2_{i+1}  (super-diagonal)
    mu  = zeros(n-1, 1);
    lam = zeros(n-1, 1);
    for i = 1 : n-1
        mu(i)  = h(i)   / (h(i) + h(i+1));
        lam(i) = h(i+1) / (h(i) + h(i+1));
    end

    % ---- branch on boundary condition type ----
    if strcmp(bc_type, 'natural')

        % natural BCs pin s2_0 = 0 and s2_n = 0 so we solve only for the
        % n-1 interior unknowns s2_1 ... s2_{n-1}.
        %
        % the (n-1) x (n-1) tridiagonal system stored as three vectors:
        %   sub-diagonal : mu(2..n-1)   length n-2
        %   main diagonal: 2 everywhere  length n-1
        %   super-diagonal: lam(1..n-2)  length n-2

        dl_nat = mu(2:n-1);          % sub-diagonal, length n-2
        d_nat  = 2 * ones(n-1, 1);  % main diagonal, length n-1
        du_nat = lam(1:n-2);        % super-diagonal, length n-2

        % solve for the interior second derivatives
        s2_int = tridiag_solve(dl_nat, d_nat, du_nat, d_int);

        % pin the endpoints to zero and glue together
        s2 = [0; s2_int; 0];

    elseif strcmp(bc_type, 'clamped')

        % clamped BCs: s'(x_0) = bc_vals(1) and s'(x_n) = bc_vals(2)
        %
        % We need two extra equations relating the endpoint first derivatives
        % to the s_i'' unknowns. Derived from differentiating the local cubic
        % on the first and last subintervals (Set 13 notes):
        %
        % Left endpoint -- differentiating p_1(x) at x_0:
        %   s'(x_0) = (f_1 - f_0)/h_1 - h_1/6 * (2*s2_0 + s2_1)
        %   rearranging:
        %   2*s2_0 + (1)*s2_1 = 6/h_1 * ((f_1-f_0)/h_1 - s'(x_0))   ... (*)
        %
        % Right endpoint -- differentiating p_n(x) at x_n:
        %   s'(x_n) = (f_n - f_{n-1})/h_n + h_n/6 * (s2_{n-1} + 2*s2_n)
        %   rearranging:
        %   (1)*s2_{n-1} + 2*s2_n = 6/h_n * (s'(x_n) - (f_n-f_{n-1})/h_n)  ... (**)
        %
        % these form the first and last rows of the full (n+1) x (n+1) system.

        d0 = (6/h(1)) * ((f_nodes(2) - f_nodes(1)) / h(1)  - bc_vals(1));
        dn = (6/h(n)) * (bc_vals(2) - (f_nodes(n+1) - f_nodes(n)) / h(n));

        % full rhs vector of length n+1
        rhs_full = [d0; d_int; dn];

        % three diagonal vectors for the (n+1) x (n+1) system:
        %   first row:     [2,  1,  0,  0, ...]  from (*)
        %   interior rows: [mu, 2, lam]           from the standard equations
        %   last row:      [..., 0, 1,  2]         from (**)

        dl_cl = [mu; 1];          % sub-diagonal length n:   mu(1..n-1) then 1
        d_cl  = 2 * ones(n+1,1); % main diagonal length n+1: all 2's
        du_cl = [1; lam];        % super-diagonal length n:  1 then lam(1..n-1)

        % solve the full system
        s2 = tridiag_solve(dl_cl, d_cl, du_cl, rhs_full);

    else
        error('bc_type must be ''natural'' or ''clamped''');
    end

    % pass knots and function values through for spline1_eval
    breaks = x_nodes;
    f_vals = f_nodes;

end