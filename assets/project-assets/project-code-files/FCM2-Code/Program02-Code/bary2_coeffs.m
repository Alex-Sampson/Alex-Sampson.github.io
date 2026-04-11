function beta = bary2_coeffs(x_nodes, n, mesh_type)
% bary2_coeffs
%   Compute barycentric Form 2 weights beta_i for the interpolating
%   polynomial using O(n) closed-form expressions that exploit the
%   structure of the node set (Berrut & Trefethen 2004, SIAM Review).
%
%   For a general mesh the weights are proportional to the Form 1 weights:
%       beta_i = gamma_i = 1 / prod_{j~=i} (x_i - x_j)
%   but for structured meshes explicit O(n) formulas are available.
%
%   The barycentric weights appear identically in both the numerator and
%   denominator of Form 2, so any common nonzero factor cancels:
%       p_n(x) = [sum_i beta_i * f_i / (x-x_i)] / [sum_i beta_i / (x-x_i)]
%   Therefore only the pattern of signs and relative magnitudes matters,
%   and the simplified (O(n)) weight formulas below are exact.
%
%   Supported mesh types and their weight formulas
%   (Berrut & Trefethen 2004):
%
%   'uniform'  -- equispaced nodes on any [a, b]  (eq. 5.1)
%       beta_i = (-1)^i * C(n, i),   i = 0, ..., n
%       where C(n,i) = n! / (i! (n-i)!) is the binomial coefficient.
%       Common factor h^n * n! is dropped since it cancels in Form 2.
%
%   'cheby1'   -- Chebyshev first kind zeros on [a, b]  (eq. 5.3)
%       x_i = cos((2i+1)*pi / (2*(n+1))),  i = 0, ..., n  (on [-1,1])
%       beta_i = (-1)^i * sin((2i+1)*pi / (2*(n+1)))
%
%   'cheby2'   -- Chebyshev second kind extrema on [a, b]  (eq. 5.4)
%       x_i = cos(i*pi/n),  i = 0, ..., n  (on [-1,1])
%       beta_i = (-1)^i * delta_i,
%       where delta_i = 1/2 for i=0 or i=n, else delta_i = 1.
%
%   'general'  -- arbitrary distinct nodes  (O(n^2) fallback)
%       beta_i = 1 / prod_{j~=i} (x_i - x_j)
%       Use this when the mesh is not one of the structured types above.
%
% Inputs:
%   x_nodes   - vector of n+1 distinct mesh points (double)
%   n         - degree of the interpolating polynomial (n+1 points)
%   mesh_type - string: 'uniform', 'cheby1', 'cheby2', or 'general'
%
% Outputs:
%   beta      - column vector of n+1 barycentric Form 2 weights (double)
%               These weights are used directly in bary2_eval in place of
%               the Form 1 gamma weights.
%
% Space:  O(n)    Time:  O(n) for structured meshes, O(n^2) for 'general'

    x_nodes = double(x_nodes(:));  % Ensure column vector in double precision
    beta    = zeros(n+1, 1);       % Pre-allocate weight vector

    signs = (-1) .^ (0:n)';        % Alternating sign vector (-1)^i, i=0..n

    switch lower(mesh_type)

        case 'uniform'
            % Berrut & Trefethen eq. 5.1: beta_i = (-1)^i * C(n,i)
            % Common scaling factors (h^n * n!) cancel in Form 2 and are dropped.
            for i = 0 : n
                beta(i+1) = signs(i+1) * nchoosek(n, i); % Binomial coefficient
            end

        case 'cheby1'
            % Berrut & Trefethen eq. 5.3: beta_i = (-1)^i * sin((2i+1)*pi/(2*(n+1)))
            % Nodes are zeros of T_{n+1}(t) on [-1,1].
            for i = 0 : n
                beta(i+1) = signs(i+1) * sin((2*i + 1) * pi / (2*(n+1)));
            end

        case 'cheby2'
            % Berrut & Trefethen eq. 5.4: beta_i = (-1)^i * delta_i
            % where delta_i = 1/2 for i=0 or i=n, else 1.
            for i = 0 : n
                if i == 0 || i == n
                    delta = 0.5;  % Half-weight for endpoint nodes
                else
                    delta = 1.0;  % Full weight for interior nodes
                end
                beta(i+1) = signs(i+1) * delta;
            end

        case 'general'
            % General O(n^2) product formula: beta_i = 1 / prod_{j~=i}(x_i - x_j)
            % Used as a fallback for arbitrary node distributions.
            for i = 1 : n+1
                prod_val = 1.0;
                for j = 1 : n+1
                    if j ~= i
                        prod_val = prod_val * (x_nodes(i) - x_nodes(j));
                    end
                end
                beta(i) = 1.0 / prod_val;
            end

        otherwise
            error('bary2_coeffs: unknown mesh_type ''%s''. Use ''uniform'', ''cheby1'', ''cheby2'', or ''general''.', mesh_type);

    end

end
