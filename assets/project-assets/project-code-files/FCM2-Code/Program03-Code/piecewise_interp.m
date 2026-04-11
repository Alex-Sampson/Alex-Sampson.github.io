function [coeffs, breaks, node_coords] = piecewise_interp(f, f_prime, a, b, M, d, node_type, hermite_flag, breaks_in)
% piecewise_interp
%   Build a piecewise polynomial g_d(x) on [a, b].
%
%   The interval [a,b] is chopped into M subintervals. By default the
%   breakpoints are uniformly spaced. Optionally a vector of M+1 custom
%   breakpoints can be passed as the 9th argument -- needed for nonuniform
%   data meshes like the one in Task 2.
%
%   If hermite_flag is true, d is forced to 3 and on each subinterval we
%   fit a cubic Hermite polynomial using f and f' at the two endpoints.
%   In that case node_type is ignored.
%
% Inputs:
%   f            - function handle for the function we're interpolating
%   f_prime      - function handle for f'(x). Only used when hermite_flag = true.
%                  Pass [] if hermite_flag = false.
%   a, b         - left and right endpoints of [a, b]
%   M            - number of subintervals
%   d            - local polynomial degree (1, 2, or 3). Overridden to 3 if hermite_flag = true.
%   node_type    - 'uniform'  : interior nodes spaced evenly inside each subinterval
%                  'cheby2'   : interior nodes placed at Chebyshev-2nd-kind points
%                  (ignored when hermite_flag = true)
%   hermite_flag - true  : piecewise cubic Hermite (needs f_prime)
%                  false : standard piecewise Lagrange (function values only)
%   breaks_in    - (optional) (M+1) x 1 vector of custom breakpoints.
%                  If omitted or empty, breakpoints are uniformly spaced.
%
% Outputs:
%   coeffs      - (d+1) x M matrix. Column s holds Newton divided-difference
%                 coefficients for the local polynomial on subinterval s.
%   breaks      - (M+1) x 1 vector of breakpoints actually used
%   node_coords - (d+1) x M matrix. Column s holds node x-coordinates
%                 for subinterval s (what newton_eval_dp needs).

    % ---- force degree to 3 when doing Hermite ----
    if hermite_flag
        d = 3;
    end

    % ---- build or accept breakpoints ----
    if nargin < 9 || isempty(breaks_in)
        % default: uniform spacing
        H = (b - a) / M;
        breaks = zeros(M+1, 1);
        for s = 0 : M
            breaks(s+1) = a + s * H;
        end
    else
        % caller supplied custom breakpoints -- use them directly
        breaks = breaks_in(:);   % ensure column vector
        if length(breaks) ~= M+1
            error('breaks_in must have exactly M+1 = %d entries', M+1);
        end
    end

    % ---- pre-allocate output matrices ----
    coeffs      = zeros(d+1, M);
    node_coords = zeros(d+1, M);

    % ---- loop over each subinterval and build the local polynomial ----
    for s = 1 : M

        a_s = breaks(s);
        b_s = breaks(s+1);

        % -- choose the d+1 interpolation nodes on [a_s, b_s] --
        if hermite_flag
            local_nodes = [a_s; b_s];   % endpoints only -- repeats handled below

        elseif strcmp(node_type, 'uniform')
            local_nodes = zeros(d+1, 1);
            for j = 0 : d
                local_nodes(j+1) = a_s + j * (b_s - a_s) / d;
            end

        elseif strcmp(node_type, 'cheby2')
            mid  = (a_s + b_s) / 2;
            half = (b_s - a_s) / 2;
            local_nodes = zeros(d+1, 1);
            for j = 0 : d
                t = cos(j * pi / d);
                local_nodes(j+1) = mid + half * t;
            end
            local_nodes = sort(local_nodes);

        else
            error('node_type must be ''uniform'' or ''cheby2''');
        end

        % -- build the Newton coefficients for this subinterval --
        if hermite_flag
            ha = a_s;
            hb = b_s;

            fa  = f(ha);
            fb  = f(hb);
            dfa = f_prime(ha);
            dfb = f_prime(hb);

            xh = [ha; ha; hb; hb];
            c  = [fa; dfa; fb; dfb];

            % order 1 divided differences
            c3_order1 = (c(3) - c(1)) / (xh(3) - xh(1));   % f[ha, hb]

            % order 2 divided differences
            c2_order2 = (c3_order1 - c(2)) / (xh(3) - xh(1));   % f[ha, ha, hb]
            c3_order2 = (c(4) - c3_order1) / (xh(4) - xh(1));   % f[ha, hb, hb]

            % order 3 divided difference
            c1_order3 = (c3_order2 - c2_order2) / (xh(4) - xh(2));   % f[ha, ha, hb, hb]

            c = [fa; dfa; c2_order2; c1_order3];

            local_nodes       = [ha; ha; hb];
            node_coords(:, s) = [ha; ha; hb; hb];

        else
            fvals = zeros(d+1, 1);
            for j = 1 : d+1
                fvals(j) = f(local_nodes(j));
            end

            c     = fvals;
            n_loc = d;

            for jj = 1 : n_loc
                for ii = n_loc : -1 : jj
                    c(ii+1) = (c(ii+1) - c(ii)) / (local_nodes(ii+1) - local_nodes(ii-jj+1));
                end
            end

            node_coords(:, s) = local_nodes;
        end

        coeffs(:, s) = c;

        if ~hermite_flag
            node_coords(:, s) = local_nodes;
        end

    end % end loop over subintervals

end