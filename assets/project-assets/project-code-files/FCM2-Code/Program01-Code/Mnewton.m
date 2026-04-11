function [x, xhist, fhist, dxhist, k] = Mnewton(f, fp, x0, m, tol, maxit)
% modified_newton
%   Solve f(x)=0 using modified Newton:
%       x_{k+1} = x_k - m*f(x_k)/f'(x_k)

    tiny = 1e-14; % Set a small number to avoid division by zero issues

    x  = x0; % Initialize x with the starting guess x0
    fx = f(x); % Evaluate the function f at the current x

    xhist  = x; % Initialize history of x values with the first guess
    fhist  = fx; % Initialize history of function values with the first evaluation
    dxhist = []; % Initialize history of changes in x as an empty array

    k = 0; % Initialize iteration counter

    while k < maxit % Loop until the maximum number of iterations is reached

        if ~isfinite(x) || ~isfinite(fx) % Check if x or fx is not a finite number
            return; % If so, exit the function
        end

        if abs(fx) < tol % Check if the absolute value of fx is less than the tolerance
            return; % If so, exit the function
        end

        dfx = fp(x); % Compute the derivative of f at the current x

        if ~isfinite(dfx) % Check if the derivative is not a finite number
            return; % If so, exit the function
        end

        if abs(dfx) < tiny % Check if the absolute value of the derivative is too small
            return; % if so, exit the function
        end

        dx = - m * fx / dfx; % Calculate the change in x using the modified Newton's method

        x  = x + dx; % Update x by adding the change
        fx = f(x); % Evaluate the function at the new x

        xhist  = [xhist; x]; % Append the new x to the history of x values
        fhist  = [fhist; fx]; % Append the new function value to the history
        dxhist = [dxhist; dx]; % Append the change in x to the history

        k = k + 1; % Increment the iteration counter

        if abs(dx) < tol*(1 + abs(x)) % Check if the change is small enough relative to x
            return; % If so, exit the function
    end
end
