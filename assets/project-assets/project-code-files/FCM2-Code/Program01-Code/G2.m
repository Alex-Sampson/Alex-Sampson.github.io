function y = G2(x)
% G2
%   Problem 6.10(c):
%       xi_{k+1}  = -sqrt(4 - eta_k^2)
%       eta_{k+1} = 1 - exp(xi_k)
%
%   Returns NaN(2,1) if the next iterate would be complex.

    xi  = x(1); % Get the first element of x, which represents xi
    eta = x(2); % Get the second element of x, which represents eta

    if 4 - eta^2 < 0 % Check if the expression under the square root is negative
        y = NaN(2,1); % If it is negative, set y to NaN to indicate an invalid result
        return; % Exit the function early since we can't proceed with a complex number
    end

    y = zeros(2,1); % Initialize y as a 2x1 vector of zeros
    y(1) = -sqrt(4 - eta^2); % Calculate the new xi value based on eta
    y(2) = 1 - exp(xi); % Calculate the new eta value based on the current xi
end
