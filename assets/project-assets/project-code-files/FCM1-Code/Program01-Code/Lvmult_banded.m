function [z] = Lvmult_banded(Lb1, Lb2, v, z, n)
% Lvmult_banded  Compute z = L*v for unit-lower L with exactly two
% lower off-diagonals stored as 1D arrays:
%   Lb1(j) = L(j+1, j),  j = 1..n-1   (first subdiagonal)
%   Lb2(j) = L(j+2, j),  j = 1..n-2   (second subdiagonal)
% Diagonal 1's are implicit.
%
% Inputs:
%   Lb1 : (n-1)-by-1
%   Lb2 : (n-2)-by-1
%   v   : n-by-1
%   z   : n-by-1 Initialized as zeros to store and overwrite
%   n   : dimension


    z = v; % Start with the diagonal contribution


  %----------------- First subdiagonal: rows j+1 get Lb1(j)*v(j)----------
  for j = 1:n-1
    z(j+1) = z(j+1) + Lb1(j) * v(j);
  end
  %-----------------------------------------------------------------------
  % --------------- Second subdiagonal: rows j+2 get Lb2(j)*v(j)----------
  for j = 1:n-2
    z(j+2) = z(j+2) + Lb2(j) * v(j);
  end
  %-----------------------------------------------------------------------
  return
end