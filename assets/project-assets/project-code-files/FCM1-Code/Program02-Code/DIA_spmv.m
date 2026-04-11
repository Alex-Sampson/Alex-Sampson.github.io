function y = DIA_spmv(DIAG, IOFF, x)
% DIA_SPMV  Compute y = C*x for a matrix C stored in DIA form.
%
% Inputs:
%   DIAG  n-by-3 matrix; column j stores the diagonal with offset IOFF(j)
%          and is padded with zeros so everything fits in one n-by-3 array
%   IOFF  1-by-3 array indicating offsets (e.g., [-k, 0, +k])
%   x     n-by-1 dense vector
%
% Output:
%   y     n-by-1 vector = C*x
%
% Rule: DIAG(i,j) multiplies x(i + IOFF(j)) when 1 <= i + IOFF(j) <= n.

    n = length(x);                  % define dimension size
    y = zeros(n,1);                 % initialize y storage
    x = x(:);                       % forces x to be a col vect.
    for j = 1:3                     % IOFF is always fixed size since we have 3 bands
        off = IOFF(j);

        if off >= 0                               % if on super diag or main diag
            imax = n - off;                       % max length of off diag array
            for i = 1:imax          
                a = DIAG(i, j);                   % define value at i,j
                if a ~= 0                         % if value is non zero
                    y(i) = y(i) + a * x(i + off); % compute the component of y
                end
            end
        else % (SUB_DIAG) off < 0: the valid rows are: i = (1-off) through n  (since i+off >= 1)
            istart = 1 - off;  
            for i = istart:n                      % *** same idea as above ****
                a = DIAG(i, j);                
                if a ~= 0
                    y(i) = y(i) + a * x(i + off);
                end
            end
        end
    end
end