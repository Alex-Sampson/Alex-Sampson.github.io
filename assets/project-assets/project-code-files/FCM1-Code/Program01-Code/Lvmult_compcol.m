function [z] = Lvmult_compcol(Lc, v, z, n)
% Lvmult_compcol  Compute w = L*v where L is unit lower triangular,
% stored in compressed-column 1D form (strict lower part).
% Lc packs column j entries for rows i=j+1..n.

  v = v(:);        % ensure column
  z = v;           % diagonal ones contribution
  k = 1;           % Initialize Index

  for j = 1:n-1
    xj = v(j);                  % Pick out entries of v at each step.
    for i = j+1:n               % multiplies the associated entry of Lc by the appropriate section of v.
      z(i) = z(i) + Lc(k) * xj; % rewrites z at each step.
      k = k + 1;
    end
  end
end