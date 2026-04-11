function [z]=Uvmult_col(M,v,z,n)

%
% Column-triad oriented version of in-place z=U*v
% with upper triangular U 
% stored in a 2d array M and vector v stored in a 1d array. M need not match dimensions with U 
% if done as a real code then M would be declared with its
% dimensions Mrdim and Mcdim (or it would come with them in an object)
% here they are extracted for reference and debugging.
% The same is true for the 1d array storing v and the 1d output array 
% storing the vector z.
% n is the dimension of U, v, and z 
% with U stored in the leading block of M in the usual fashion.
%
% This version would be preferred in a language that
% stores two dimensional arrays in a column major ordering.
% In a language that uses row major ordering for 
% two dimensional arrays, a row oriented inner product
% form of algorithm would be preferred.
%
% The output vector z  would be declared in a real code.
% It is passed here as an input variable to indicate
% how many languages would provide a data structure 
% to be updated with the computed output, e.g.,
% the product U*v could go straight into z or 
% you might update it z <-- z + U*v
% Here its dimensions are not needed but they would be passed in
% and used in that declaration.
% Here they are extracted for reference and debugging.
%
%
%
% The following extract the dimensions
% of the M, v and z arrays from the information
% Matlab maintains for them in their particular
% objects. They would be used in the declarations
% of M, v and z mentioned above in a code with 
% a compiled and typed language.
% 

  Mrdim=size(M,1);
  Mcdim=size(M,2); 
  zdim=size(z,1);
  vdim=size(v,1);  
%
% note the code relies on Matlab suppressing operations on 
% vectors of length 0, i.e., when index range r:s is such that r > s
%


%
% the i loop goes over columns of U to scale them with an element of v
% and contribute to a running sum.
%
% Column i of U has its possibly nonzero elements in positions  (1:i,i).
% all positions (i,j) with i > j are 0 and not stored.
%


  z=zeros(n,1);  % put zeros only in the elements corresponding to the
                 % mathematical dimension not the entire array
  for i = 1:n
    z(1:i) = z(1:i) + M(1:i,i)*v(i);
  end % end i loop
  return