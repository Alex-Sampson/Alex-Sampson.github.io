function y = f3(x, x_nodes)
% f3
%   Evaluate f3(x) = l_n(x), the n-th Lagrange basis function for the
%   mesh x_nodes, in single precision.
%
%       f3(x_i) = 0,  0 <= i <= n-1
%       f3(x_n) = 1
%
%   Computed in product form:
%
%       l_n(x) = prod_{i=0}^{n-1} (x - x_nodes(i))
%                / prod_{i=0}^{n-1} (x_nodes(n) - x_nodes(i))
%
% Inputs:
%   x        - vector of query points (single)
%   x_nodes  - vector of n+1 mesh points; l_n is the basis function
%              associated with the last node x_nodes(n+1) (single)
%
% Outputs:
%   y        - vector of function values l_n(x), same size as x (single)
%
% Space:  O(length(x))    Time:  O(n * length(x))

    x       = single(x);        % Ensure query points are in single precision
    x_nodes = single(x_nodes);  % Ensure mesh points are in single precision

    n   = length(x_nodes) - 1;  % Degree of the polynomial
    x_n = x_nodes(n+1);         % The node associated with this basis function

    numer = ones(size(x), 'single'); % Initialize numerator product
    denom = single(1.0);             % Initialize denominator scalar

    for i = 1 : n                            % Loop over all nodes except the last
        numer = numer .* (x - x_nodes(i));   % Accumulate (x - x_i) into numerator
        denom = denom  * (x_n - x_nodes(i)); % Accumulate (x_n - x_i) into denominator
    end

    y = numer / denom; % Divide to get the Lagrange basis function values

end
