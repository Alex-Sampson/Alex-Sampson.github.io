function x_leja = leja_order(x)
% leja_order
%   Reorder mesh points into Leja ordering.
%
%       x_leja(1)   = argmax_i |x(i)|
%       x_leja(k+1) = argmax_{x not yet chosen} prod_{j=1}^{k} |x - x_leja(j)|
%
% Inputs:
%   x       - vector of n+1 distinct mesh points
%
% Outputs:
%   x_leja  - column vector of the same n+1 points in Leja order
%
% Space:  O(n)    Time:  O(n^2)

    x        = x(:);            % Force x to a column vector
    m        = length(x);       % Total number of points
    x_leja   = zeros(m, 1);     % Pre-allocate output vector
    selected = false(m, 1);     % Track which points have been selected

    [~, idx]      = max(abs(x));  % Find index of largest magnitude point
    x_leja(1)     = x(idx);       % Set first Leja point
    selected(idx) = true;          % Mark it as selected

    for k = 1 : m-1                        % Loop to select each remaining point
        node_prod            = ones(m, 1); % Initialize running node product for all candidates
        for j = 1 : k                      % Accumulate product over already-selected points
            node_prod = node_prod .* abs(x - x_leja(j)); % Multiply in |x - x_leja(j)|
        end
        node_prod(selected)   = -1;          % Mask selected points from consideration
        [~, idx]              = max(node_prod); % Find candidate maximising node product
        x_leja(k+1)           = x(idx);        % Store next Leja point
        selected(idx)         = true;           % Mark it as selected
    end

end
