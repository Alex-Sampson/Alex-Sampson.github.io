function [x, xhist, fhist, dxhist, k, a, b] = regula_falsi(f, a0, b0, tol, maxit)
% regula_falsi
%   Solve f(x)=0 using false position (Regula Falsi).
%   Maintains a bracket [a,b] with f(a)*f(b) < 0.
%
% Outputs:
%   x      - final iterate (or NaN on failure)
%   xhist  - iterate history
%   fhist  - f(x) history
%   dxhist - step history (x_k - x_{k-1})
%   k      - number of iterations actually taken (matches length(xhist))
%   a,b    - final bracket endpoints

    tiny = 1e-14; % A tollerance to avoid division by zero issues


    x      = NaN; % Initialize x to NaN, will hold the final result
    xhist  = [];  % Initialize history of x values
    fhist  = [];  % Initialize history of f(x) values
    dxhist = [];  % Initialize history of step sizes
    k      = 0;   % Initialize iteration counter

    a = a0; % Set the left endpoint of the bracket
    b = b0; % Set the right endpoint of the bracket


    if ~isfinite(a) || ~isfinite(b) % Check if a and b are finite numbers
        return; % If not, exit the function
    end

    fa = f(a); % Evaluate the function at the left endpoint
    fb = f(b); % Evaluate the function at the right endpoint

    if ~isfinite(fa) || ~isfinite(fb) % Check if function values are finite
        return; % If not, exit the function
    end

   
    if abs(fa) < tol % Check if the left endpoint is a root
        x = a; xhist = a; fhist = fa; k = 1; % Set outputs accordingly
        return; % Exit the function
    end
    if abs(fb) < tol % Check if the right endpoint is a root
        x = b; xhist = b; fhist = fb; k = 1; % Set outputs accordingly
        return; % Exit the function
    end

  
    if fa*fb > 0 % Check if the function values have the same sign
        return; % If they do, exit the function as we don't have a valid bracket
    end

    xprev = NaN; % Initialize previous x value to NaN

    while k < maxit % Loop until we reach the maximum number of iterations

        denom = fb - fa; % Calculate the denominator for the false position formula
        if ~isfinite(denom) || abs(denom) < tiny % Check for stability
            return; % If not stable, exit the function
        end

       
        x  = (a*fb - b*fa) / denom; % Calculate the new x using the false position formula
        fx = f(x); % Evaluate the function at the new x

        if ~isfinite(x) || ~isfinite(fx) % Check if x and f(x) are finite
            return; % If not, exit the function
        end

      
        xhist = [xhist; x]; % Append the new x to the history
        fhist = [fhist; fx]; % Append the new f(x) to the history
        if ~isnan(xprev) % Check if we have a previous x value
            dxhist = [dxhist; x - xprev]; % Calculate and store the step size
        end

        
        k = k + 1; % Increment the iteration counter

      
        if abs(fx) < tol % Check if we found a root
            return; % If yes, exit the function
        end
        if ~isnan(xprev) % Check if we have a previous x value
            if abs(x - xprev) < tol*(1 + abs(x)) % Check for convergence
                return; % If converged, exit the function
            end
        end

       
        if fa*fx < 0 % Check if the root is in the left half
            b  = x; % Update the right endpoint
            fb = fx; % Update the function value at the right endpoint
        else % Otherwise, the root is in the right half
            a  = x; % Update the left endpoint
            fa = fx; % Update the function value at the left endpoint
        end

        xprev = x; % Update the previous x value for the next iteration
    end
end
