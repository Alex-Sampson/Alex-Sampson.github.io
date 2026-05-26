function [Q, n_evals] = comp_gauss2(f, a, b, M)
% comp_gauss2
%   Composite 2-point Gauss-Legendre Rule on [a, b] with M uniform subintervals.
%   This function breaks the interval [a,b] into M equal pieces and on each
%   piece uses two smart sample points (Gauss nodes) to approximate the
%   integral. The rule is accurate for cubic polynomials and usually very
%   good for smooth functions. It always evaluates the function twice per
%   subinterval, so total function calls = 2*M.
%
% Inputs:
%   f    - function handle @(x) ... representing the integrand
%   a, b - left and right endpoints of the integration interval
%   M    - number of equal subintervals to use (positive integer)
%
% Outputs:
%   Q       - the approximate integral over [a,b]
%   n_evals - number of function evaluations performed (2*M)

    H = (b - a) / M;   % compute the width of each small subinterval

    % -- 2-point Gauss-Legendre nodes and weights on [-1, 1] --
    % The two sample points on the standard interval are +/- 1/sqrt(3),
    % and both carry weight 1. These are fixed numbers for the 2-point rule.
    t1 = -1.0 / sqrt(3);   % left reference node on [-1,1]
    t2 =  1.0 / sqrt(3);   % right reference node on [-1,1]

    % -- loop over each subinterval and accumulate --
    Q = 0;                 % initialize the running total for the integral
    for k = 1 : M          % process each subinterval from 1 to M

        % midpoint of subinterval k: a + (k - 1/2) * H
        m_k = a + (k - 0.5) * H;   % center point of the current subinterval

        % half-width used in the affine map [-1,1] -> [a_k, b_k]
        half_H = H / 2;            % half the subinterval length (Jacobian factor)

        % map the two reference nodes into the physical subinterval
        x1 = m_k + t1 * half_H;    % first physical Gauss node in subinterval k
        x2 = m_k + t2 * half_H;    % second physical Gauss node in subinterval k

        % local contribution: Jacobian * (w1*f(x1) + w2*f(x2))
        % weights are both 1 so this simplifies to half_H * (f(x1) + f(x2))
        Q = Q + half_H * (f(x1) + f(x2));  % add the two-sample estimate for this piece

    end                     % end loop over subintervals

    % -- exactly 2 evaluations per subinterval, no sharing --
    n_evals = 2 * M;        % record how many function calls were made

end