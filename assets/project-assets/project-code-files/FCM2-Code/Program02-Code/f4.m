function y = f4(x)
% f4
%   Evaluate f4(x) = 1 / (1 + 25*x^2) (single precision).
%
%   This is the Runge function, commonly used to demonstrate the
%   failure of high-degree polynomial interpolation on uniform meshes.
%   No product form is available for this function.
%
% Inputs:
%   x    - vector of query points (single)
%
% Outputs:
%   y    - vector of function values 1/(1+25x^2), same size as x (single)
%
% Space:  O(length(x))    Time:  O(length(x))

    x = single(x);                          % Ensure query points are in single precision
    y = single(1.0) ./ (single(1.0) + single(25.0) .* x .^ 2); % Evaluate Runge function

end
