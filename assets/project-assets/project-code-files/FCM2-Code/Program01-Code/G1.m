function y = G1(x)
% G1
%   Problem 6.10(b):
%       xi_{k+1}  = log(1 - eta_k)
%       eta_{k+1} = -sqrt(4 - xi_k^2)
%
%   Returns NaN(2,1) if the next iterate would be complex.

    xi  = x(1); % Get the first element of x and assign it to xi
    eta = x(2); % Get the second element of x and assign it to eta

    if 1 - eta <= 0 % Check if 1 minus eta is less than or equal to zero
        y = NaN(2,1); % If true, set y to a 2x1 array of NaNs
        return; % Exit the function early since we can't proceed
    end

    if 4 - xi^2 < 0 % Check if 4 minus the square of xi is less than zero
        y = NaN(2,1); % If true, set y to a 2x1 array of NaNs
        return; % Exit the function early since we can't proceed
    end

    y = zeros(2,1); % Initialize y as a 2x1 array of zeros
    y(1) = log(1 - eta); % Calculate the logarithm of (1 - eta) and store it in the first element of y
    y(2) = -sqrt(4 - xi^2); % Calculate the negative square root of (4 - xi^2) and store it in the second element of y
end
