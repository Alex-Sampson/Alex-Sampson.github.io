function [Acomp] = LUmult_outer(M, n, mode)

% M - Packed matrix storing L and U
% n - size of matrix
% mode - if mode=1 we are doing |L||U|, else LU.
%
% Compute Acomp = L*U from packed storage M, without unpacking.
% Packing in M:
%   - i > j : M(i,j) = L(i,j)          (strict lower of L)
%   - i <= j: M(i,j) = U(i,j)          (diag + upper of U)
%   - L has unit diagonal: L(i,i) = 1  (implicit, not stored)
%
% Outer-product scheme:
%   For k = 1..n, add L(:,k)*U(k,:) but only where entries can be nonzero:
%     rows i = k..n  (since L(i,k)=0 for i<k; L(k,k)=1)
%     cols j = k..n  (since U(k,j)=0 for j<k)

if mode==1
    M = abs(M);
end 

  Acomp=zeros(n); %Initializes A as an n-by-nx of zeros.

  for k = 1:n     % iterate rows where L(i,k) may be nonzero
    for i = k:n   % L(i,k): 1 on the diagonal, otherwise from strict-lower of M
      if i == k   % diag elements of L
        l = 1;    
      else
        l = M(i,k);      % L(i,k) from M
      end

      if l ~= 0          % iterate cols where U(k,j) may be nonzero
        for j = k:n      % U(k,j) sits in the upper/diag part of M
          u = M(k,j);    % U(k,j) from packed M
          if u ~= 0
            Acomp(i,j) = Acomp(i,j) + l * u;
          end
        end
      end
    end
  end

  return
end




