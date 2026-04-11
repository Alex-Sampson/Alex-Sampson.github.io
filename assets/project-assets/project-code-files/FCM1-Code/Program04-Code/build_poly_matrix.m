function B = build_poly_matrix(x, d, basis)
% BUILD_POLY_MATRIX  Build a polynomial regression matrix B for a chosen basis.
%
% Given x = (x_1, ..., x_n)', and a nonnegative integer d, this forms
% the n-by-(d+1) matrix
%
%       B = [phi_0(x), phi_1(x), ..., phi_d(x)],
%
% where {phi_k} is either the monomial basis or the Chebyshev T_k basis.
%
% INPUT:
%   x     : n-by-1 vector of evaluation points (if row, we convert to column)
%   d     : nonnegative integer, maximum degree
%   basis : integer flag:
%             0---> monomial basis
%             1---> Chebyshev basis (first kind)
%
% OUTPUT:
%   B     : n-by-(d+1) matrix. Column j stores phi_{j-1}(x) for j = 1,...,d+1.
% ===================================================================================================

   
    x = x(:); % Make sure x is a column vector
    n = length(x);

%--------------------------------Error Message-------------------------------------------------------

    if d < 0 || floor(d) ~= d
        error('build_poly_matrix: degree d must be a nonnegative integer.');
    end

    if basis ~= 0 && basis ~= 1
        error('build_poly_matrix: basis must be 0 (monomial) or 1 (Chebyshev).');
    end
    
%----------------------------------------------------------------------------------------------------


    % Initialize B
    B = zeros(n, d+1);

%------------------------------ MONOMIAL BASIS (basis = 0) ------------------------------------------

    if basis == 0

        % Column 1: phi_0(x) = 1
        for i = 1:n
            B(i, 1) = 1;
        end

        
 % ---------------- For each row i, build 1, x(i), x(i)^2, ..., x(i)^d ------------------------------
 
        for i = 1:n
            p = 1;                 % will hold x(i)^k
            for k = 1:d
                p = p * x(i);      % update p from x(i)^(k-1) to x(i)^k
                B(i, k+1) = p;     % store in column k+1
            end
        end

% ----------------------------- CHEBYSHEV BASIS (basis = 1) -----------------------------------------

    else  % basis == 1

        % Chebyshev T_k basis:
        %   T_0(x) = 1
        %   T_1(x) = x
        %   T_{m+1}(x) = 2*x*T_m(x) - T_{m-1}(x)
        %
        % We store:
        %   column 1 = T_0
        %   column 2 = T_1
        %   column 3 = T_2
        %   ...
        %   column d+1 = T_d

        % Column 1: T_0(x) = 1
        for i = 1:n
            B(i, 1) = 1;
        end

        % Column 2: T_1(x) = x (only if d >= 1)
        if d >= 1
            for i = 1:n
                B(i, 2) = x(i);
            end
        end

        % Use the recurrence for m = 1,...,d-1 to build T_{m+1}
        %   T_{m+1}(x_i) = 2 * x_i * T_m(x_i) - T_{m-1}(x_i)
        %
        % Here:
        %   column m+1 stores T_m,
        %   column m   stores T_{m-1},
        %   we write T_{m+1} into column (m+2).
        
        for m = 1:(d-1)
            for i = 1:n
                Tm   = B(i, m+1);                  % T_m(x_i)
                Tm_1 = B(i, m);                    % T_{m-1}(x_i)
                B(i, m+2) = 2 * x(i) * Tm - Tm_1;
            end
        end

    end
end