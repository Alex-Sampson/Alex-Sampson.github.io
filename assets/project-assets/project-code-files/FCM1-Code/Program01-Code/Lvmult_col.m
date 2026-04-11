function [z]=Lvmult_col(M,v,z,n)

%A function that takes in arguments:
% M -- The UNIT Lower triangular matrix doing the multiplication
% v -- The vector that we are operating on
% z -- a vector initialized as v, that we will update at each step
% of a loop -- again, we technically dont need this as an argument, its
% defined inside the function, it will just stay for driver compatability
% purposes.

% n -- the dimension of our matrix/vector's (M is square since it's lower
% triangular)



  Mrdim=size(M,1); % Tells us how many rows M has
  Mcdim=size(M,2); % Tells us how many collumns M has 
                   % If M is square, Mcdim=Mrdim
  zdim=size(z,1);  % Tells us the length of z, note z is stored as a row vector
  vdim=size(v,1);  % Tells us the length of v, note v is stored as a row vector



  z = v(1:n);                                 %Initializes z as v
for j = 1:n-1
  z(j+1:n) = z(j+1:n) + M(j+1:n, j) * v(j);   %updates the 'tails' of z at each step by
end                                           %adding the next col mult. of L with v