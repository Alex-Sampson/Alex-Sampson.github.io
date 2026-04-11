function H_n = hilbert_H(kappa_y)
% hilbert_H
%   Compute the H_n norm summary statistic from the condition number
%   kappa(x, n, y) evaluated at a vector of query points.
%
%       H_n = max_{a <= x <= b} kappa(x, n, y)
%
%   H_n summarises the worst-case amplification of data errors in the
%   function values y_i when evaluating the interpolating polynomial.
%   It is approximated as the maximum of kappa(x,n,y) over the supplied
%   query points.  The accuracy of this approximation depends on the
%   density of the query grid passed to cond_kappa.
%
% Inputs:
%   kappa_y  - vector of kappa(x,n,y) values at query points,
%              as returned by cond_kappa (double)
%
% Outputs:
%   H_n      - scalar approximation of the H_n norm (double)
%
% Space:  O(1)    Time:  O(length(kappa_y))

    H_n = max(kappa_y); % H_n is the maximum of kappa(x,n,y) over all query points

end
