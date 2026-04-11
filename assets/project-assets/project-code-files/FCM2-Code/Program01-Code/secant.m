function [x, xhist, fhist, dxhist, k] = secant(f, x0, x1, tol, maxit)
% secant
%   Solve f(x)=0 using the secant method:
%       x_{k+1} = x_k - f(x_k)*(x_k-x_{k-1})/(f(x_k)-f(x_{k-1}))

    tiny = 1e-14; % A tollerance to avoid division by zero

    xm1 = x0; % Previous guess
    x   = x1; % Current guess

    fm1 = f(xm1); % Function value at the previous guess
    fx  = f(x); % Function value at the current guess

    xhist  = [xm1; x]; % History of x values
    fhist  = [fm1; fx]; % History of function values
    dxhist = [x - xm1]; % Change in x values

    k = 0; % Iteration counter

    while k < maxit % Loop until we reach the maximum number of iterations

        if ~isfinite(x) || ~isfinite(fx) || ~isfinite(xm1) || ~isfinite(fm1)
            return; % Exit if any value is not finite (like Inf or NaN)
        end

        if abs(fx) < tol
            return; % Exit if the function value is within the tolerance
        end

        denom = fx - fm1; % Calculate the denominator for the secant formula

        if ~isfinite(denom)
            return; % Exit if the denominator is not finite
        end

        if abs(denom) < tiny
            return; % Exit if the denominator is too small (to avoid division by zero)
        end

        xnew = x - fx*(x - xm1)/denom; % Calculate the new guess using the secant formula
        dx   = xnew - x; % Change in x

        xm1 = x; % Update previous guess to current guess
        fm1 = fx; % Update previous function value to current function value

        x  = xnew; % Move to the new guess
        fx = f(x); % Calculate the function value at the new guess

        xhist  = [xhist; x]; % Append the new guess to the history
        fhist  = [fhist; fx]; % Append the new function value to the history
        dxhist = [dxhist; dx]; % Append the change in x to the history

        k = k + 1; % Increment the iteration counter

        if abs(dx) < tol*(1 + abs(x))
            return; % Exit if the change in x is small enough
        end
    end
end
