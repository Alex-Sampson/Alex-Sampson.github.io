function [v, tau] = hh_vec(x)
% HH_VEC  Builds the vector and scalar needed in the construction of each householder reflector using each column x of A.
%
%
%   Returns v and tau so that we can begin to construct the matrix

%       H = I - tau * v * v'

%   Such that:

%       H * x = [alpha; 0; 0; ...],   alpha = norm(x,2) >= 0.

%   Note: v = 1/u_1 (u) where u = x-(alpha)e_1 where e_1 is the basis vector in R^{n-k} at step k we            
%         scale it so we dont have to store it explicitly
%
% INPUT:
%   x : The column vector of A at the i-th transformation step
%
% OUTPUT:
%   v   : Householder vector of the same length as x, with v(1) = 1.
%   tau : nonnegative scalar so that H = I - tau * v * v' is the reflector.
%
% -------------------------------------------------------------------------

    % Make sure x is a COLUMN vector 
    x = x(:);
    m = length(x);

    % Initialize v and tau, these will also be default outputs (identity reflector) in case x is zero
    v   = zeros(m,1);
    tau = 0;

    % alpha = ||x||_2  
    alpha = norm(x, 2);

    if alpha == 0 % x is the zero vector: H = I works (do nothing)
        v(1) = 1;
        return
    end

%-----------------------------------------------------------------------------------------------------
    % We want H*x = [alpha; 0; ...] with alpha >= 0.
    % Form u = x - alpha*e1, then scale so v(1)=1. I.e. v=1/(u_1)*u.
%-----------------------------------------------------------------------------------------------------

    u    = x;
    u(1) = u(1) - alpha;   % careful: could be 0 if x = [alpha; 0; ...].

%-----------------------------------------------------------------------------------------------------
%    If u(1)=0 then v=1/(u_1)*u is undefined, instead since u=(\alpha,0,...,0) we can just use the 
%    identity transformation.
%-----------------------------------------------------------------------------------------------------
    if u(1) == 0
        v(1) = 1;
        tau  = 0;
        return
    end

%-----------------------------------------------------------------------------------------------------
%    Otherwise compute v as specified above
%-----------------------------------------------------------------------------------------------------

    v = u / u(1);

%-----------------------------------------------------------------------------------------------------
%    Compute tau
%-----------------------------------------------------------------------------------------------------

    vv  = norm(v,2)^2;
    tau = 2 / vv;
    
    
end