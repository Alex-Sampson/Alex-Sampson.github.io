function [x, xhist, fhist, dxhist, k] = steffensen(f, x0, tol, maxit)
% steffensen
%   Solve f(x)=0 using Steffensen's method (derivative-free):
%       x_{k+1} = x_k - f(x_k)^2 / ( f(x_k+f(x_k)) - f(x_k) )

    tiny = 1e-14; % A tollerance to prevent division by zero

    x  = x0; % Initialize x with the starting guess
    fx = f(x); % Evaluate the function at the starting guess

    xhist  = x; % Initialize history of x values with the starting guess
    fhist  = fx; % Initialize history of function values with the initial function value
    dxhist = []; % Initialize history of changes in x

    k = 0; % Initialize iteration counter

    while k < maxit % Loop until the maximum number of iterations is reached

        if ~isfinite(x) || ~isfinite(fx) % Check if x or fx is not a finite number
            return; % Exit if we hit a non-finite value
        end

        if abs(fx) < tol % Check if the function value is within the tolerance
            return; % Exit if we have a good enough solution
        end

        y  = x + fx; % Calculate the next point y based on current x and fx
        fy = f(y); % Evaluate the function at the new point y

        if ~isfinite(y) || ~isfinite(fy) % Check if y or fy is not a finite number
            return; % Exit if we hit a non-finite value
        end

        denom = fy - fx; % Calculate the denominator for the update formula

        if ~isfinite(denom) % Check if the denominator is not a finite number
            return; % Exit if we hit a non-finite value
        end

        if abs(denom) < tiny % Check if the denominator is too small
            return; % Exit to avoid division by zero
        end

        xnew = x - fx*fx/denom; % Update x using Steffensen's formula
        dx   = xnew - x; % Calculate the change in x

        x  = xnew; % Update x to the new value
        fx = f(x); % Evaluate the function at the new x

        xhist  = [xhist; x]; % Append the new x to the history
        fhist  = [fhist; fx]; % Append the new function value to the history
        dxhist = [dxhist; dx]; % Append the change in x to the history

        k = k + 1; % Increment the iteration counter

        if abs(dx) < tol*(1 + abs(x)) % Check if the change is within the tolerance
            return; % Exit if we have converged
        end
    end
end
