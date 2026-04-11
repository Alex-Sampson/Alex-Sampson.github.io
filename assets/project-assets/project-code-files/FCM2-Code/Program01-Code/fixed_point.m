function [x, xhist, rhist, dxhist, k] = fixed_point(phi, x0, tol, maxit)
% fixed_point
%   Fixed-point iteration:
%       x_{k+1} = phi(x_k)
%   Records r_k = x_k - phi(x_k)

    x = x0; % Start with the initial guess x0

    px = phi(x); % Apply the function phi to the current guess x
    r  = x - px; % Calculate the residual, which is the difference between x and phi(x)

    xhist  = x; % Initialize the history of x values with the first guess
    rhist  = r; % Initialize the history of residuals with the first residual
    dxhist = []; % Initialize the history of changes in x as an empty array

    k = 0; % Start the iteration count at 0

    while k < maxit % Continue iterating until we reach the maximum number of iterations

        if ~isfinite(x) || ~isfinite(px) || ~isfinite(r) % Check if any of the values are not finite
            return; % If any value is not finite, exit the function
        end

        if abs(r) < tol % Check if the absolute value of the residual is less than the tolerance
            return; % If it is, we can stop iterating as we've converged
        end

        xnew = px; % Update the new guess to the current value of phi(x)
        dx   = xnew - x; % Calculate the change in x

        x  = xnew; % Set x to the new guess
        px = phi(x); % Calculate the new value of phi at the updated x
        r  = x - px; % Update the residual with the new values

        xhist  = [xhist; x]; % Append the new x value to the history
        rhist  = [rhist; r]; % Append the new residual to the history
        dxhist = [dxhist; dx]; % Append the change in x to the history

        k = k + 1; % Increment the iteration count

        if abs(dx) < tol*(1 + abs(x)) % Check if the change in x is small enough relative to the tolerance
            return; % If it is, we can stop iterating as we've converged
        end
    end
end
