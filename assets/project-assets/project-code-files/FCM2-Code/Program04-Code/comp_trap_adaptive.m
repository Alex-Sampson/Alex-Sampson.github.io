function [Q, M_final, n_evals] = comp_trap_adaptive(f, a, b, tol, M0)
% comp_trap_adaptive
%   Composite Trapezoidal Rule with adaptive step-halving (alpha = 1/2).
%   Starts with M0 subintervals and doubles the number of subintervals
%   until a simple error estimate says the result is within tol.
%
%   When we double the number of intervals, the old sample points are
%   still valid. We only need to evaluate f at the new midpoints that
%   appear between the old nodes. That saves work.
%
%   The recurrence used to update the trapezoid result when doubling
%   the intervals is:
%       I_{2n} = (1/2) * I_n + H_next * sum(f at new midpoints)
%   where H_next is the new step size.
%
%   We estimate the error using the fact that the trapezoid rule's
%   main error term scales like H^2. Halving H makes that term about
%   1/4 as big, so the difference between the fine and coarse
%   results gives a cheap error estimate:
%       estimated_error ≈ (1/3) * |I_fine - I_coarse|
%
% Inputs:
%   f    - function handle @(x) ...
%   a, b - left and right endpoints
%   tol  - error tolerance -- stop when estimated error < tol
%   M0   - initial number of subintervals (optional, default = 1)
%
% Outputs:
%   Q       - best approximation when stopping criterion was met
%   M_final - number of subintervals used in the final approximation
%   n_evals - total function evaluations across ALL refinement levels

    % -- handle the optional M0 argument --
    if nargin < 5                       % if the caller didn't give M0
        M0 = 1;                         % use 1 subinterval to start
    end                                  % end the nargin check

    % -- compute the initial coarse approximation using comp_trap --
    % this gives us I_coarse and burns M0+1 function evaluations
    [Q_coarse, evals_init] = comp_trap(f, a, b, M0); % call comp_trap once
    n_evals = evals_init;                % record the evaluations used so far
    M_curr  = M0;                        % set current number of intervals

    % -- refinement loop: keep halving H until estimated error < tol --
    while true                           % loop until we break explicitly

        % the next level has twice as many subintervals
        M_next = 2 * M_curr;             % new number of intervals = 2 * current
        H_next = (b - a) / M_next;       % compute the new (halved) step size

        % the new midpoints are the centers of each current subinterval.
        % for subinterval k (1-indexed), the midpoint is:
        %   a + (2k - 1) * H_next   for k = 1, 2, ..., M_curr
        % these are the ONLY points we need to evaluate f at this level.
        midpoint_sum = 0;                % initialize accumulator for new f-values
        for k = 1 : M_curr               % loop over the number of new midpoints
            x_mid = a + (2*k - 1) * H_next; % compute the k-th midpoint location
            midpoint_sum = midpoint_sum + f(x_mid); % add f at this midpoint
        end                              % end midpoint loop

        % M_curr new evaluations at this level -- that's it
        n_evals = n_evals + M_curr;      % update total function evaluation count

        % -- apply the recycling recurrence from Set 22 --
        % I_{2n} = (1/2)*I_n + H_next * sum_of_new_midpoints
        % this is algebraically identical to running comp_trap fresh with
        % M_next intervals, but avoids re-evaluating the existing nodes
        Q_fine = Q_coarse / 2 + H_next * midpoint_sum; % compute refined integral

        % -- Richardson error estimate --
        % for O(H^2) methods the leading error scales as 1/4 when H halves,
        % so (Q_fine - Q_coarse) ~ (3/4) * true_error_fine, rearranging:
        %   true_error_fine ~ (1/3) * |Q_fine - Q_coarse|
        err_est = abs(Q_fine - Q_coarse) / 3; % estimate the true error on fine grid

        % -- update state for the next possible iteration --
        Q_coarse = Q_fine;               % promote fine result to current coarse
        M_curr   = M_next;               % update current number of intervals

        % -- check stopping criterion --
        if err_est < tol                 % if estimated error below tolerance
            break                        % we are done; exit the loop
        end                              % end tolerance check

        % safety valve -- if something is going wrong and M explodes, bail out
        if M_curr > 2^20                 % if intervals exceed a very large number
            warning('comp_trap_adaptive: hit M = 2^20 without converging. Returning current estimate.'); % warn user
            break                        % exit the loop to avoid runaway compute
        end                              % end safety check

    end                                  % end adaptive refinement loop

    Q       = Q_coarse;                   % final approximation to return
    M_final = M_curr;                      % final number of subintervals used

end