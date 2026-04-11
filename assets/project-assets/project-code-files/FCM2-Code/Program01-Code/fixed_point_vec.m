function [x, xhist, rhist, dxhist, k] = fixed_point_vec(G, x0, tol, maxit)
% fixed_point_vec
%   Vector fixed-point iteration:
%       x_{k+1} = G(x_k)
%   Records:
%       r_k  = ||x_k - G(x_k)||_2
%       dx_k = ||x_{k+1} - x_k||_2
%
%   Stops if G(x) produces NaN or Inf.

    x = x0; % Start with the initial guess x0
    Gx = G(x); % Compute G(x) for the initial guess

    if any(~isfinite(Gx)) % Check if G(x) is not finite (NaN or Inf)
        x = NaN(size(x0)); % If so, set x to NaN of the same size as x0
        xhist = x0; % Record the initial guess in the history
        rhist = NaN; % Set the residual history to NaN
        dxhist = []; % Initialize the change history as empty
        k = 0; % Set the iteration count to 0
        return; % Exit the function early
    end

    r = norm(x - Gx); % Calculate the norm of the difference between x and G(x)

    xhist  = x; % Initialize the history of x with the initial guess
    rhist  = r; % Store the initial residual in the history
    dxhist = []; % Initialize the change history as empty

    k = 0; % Start the iteration count at 0

    while k < maxit % Loop until we reach the maximum number of iterations

        if any(~isfinite(x)) || ~isfinite(r) % Check if current x or residual is not finite
            return; % If so, exit the function
        end

        if r < tol % Check if the residual is less than the tolerance
            return; % If it is, we can stop the iteration
        end

        xnew = Gx; % Update xnew to the current G(x)
        dx   = norm(xnew - x); % Calculate the change between the new and old x

        x  = xnew; % Update x to the new value
        Gx = G(x); % Compute G(x) for the updated x

        if any(~isfinite(Gx)) % Check if the new G(x) is not finite
            return; % If so, exit the function
        end

        r = norm(x - Gx); % Calculate the new residual

        xhist  = [xhist, x]; % Append the new x to the history
        rhist  = [rhist; r]; % Append the new residual to the history
        dxhist = [dxhist; dx]; % Append the change to the history

        k = k + 1; % Increment the iteration count

        if dx < tol*(1 + norm(x)) % Check if the change is less than the tolerance
            return; % If it is, we can stop the iteration
        end
    end
end
