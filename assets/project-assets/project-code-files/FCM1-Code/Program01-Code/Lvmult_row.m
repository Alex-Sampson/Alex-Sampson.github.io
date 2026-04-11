function [z]=Lvmult_row(M,v,z,n)

%A function that takes in arguments:
% M -- The Lower triangular matrix doing the multiplication
% v -- The vector that we are operating on
% z -- a vector initialized as all zeros that we will update at each step
% of a loop -- this is not used for this code, we could drop z, from the
% function variables, but for uniformity we will keep it. INFACT we DONT
% need it AT ALL, I think we just keep if for Driver compatability.

% n -- the dimension of our matrix/vector's (M is square since it's lower
% triangular)

  Mrdim=size(M,1); % Tells us how many rows M has
  Mcdim=size(M,2); % Tells us how many collumns M has 
                   % If M is square, Mcdim=Mrdim
  zdim=size(z,1);  % Tells us the length of z, note z is stored as a row vector
  vdim=size(v,1);  % Tells us the length of v, note v is stored as a row vector

  
for i = 1:n
  z(i) = v(i) + M(i,1:i-1) * v(1:i-1);   % when i=1, the product is [], returns 0
end
  
