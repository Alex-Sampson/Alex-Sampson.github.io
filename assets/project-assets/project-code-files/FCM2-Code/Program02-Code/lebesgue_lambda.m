function Lambda_n = lebesgue_lambda(kappa_1)
% lebesgue_lambda
%   Compute the Lebesgue constant Lambda_n from the Lebesgue function
%   evaluated at a vector of query points.
%
%       Lambda_n = max_{a <= x <= b} kappa(x, n, 1)
%                = max_{a <= x <= b} sum_{i=0}^{n} |l_i(x)|
%
%   The Lebesgue constant is approximated as the maximum of kappa(x,n,1)
%   over the supplied query points.  The accuracy of this approximation
%   depends on the density of the query grid passed to cond_kappa.
%
% Inputs:
%   kappa_1   - vector of Lebesgue function values at query points,
%               as returned by cond_kappa (double)
%
% Outputs:
%   Lambda_n  - scalar approximation of the Lebesgue constant (double)
%
% Space:  O(1)    Time:  O(length(kappa_1))

    Lambda_n = max(kappa_1); % Lebesgue constant is the maximum of the Lebesgue function

end
