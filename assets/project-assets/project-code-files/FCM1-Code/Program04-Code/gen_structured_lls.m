function [A, b, xtrue, Rtrue, ctrue, dtrue] = gen_structured_lls(n, k, kind, dscale, nref)
    
% GEN_STRUCTURED_LLS  Generate structured least-squares test problems.
%
% We first build a "nice" problem in an orthogonal coordinate system:
%
%        [Rtrue] x ≈ [ctrue]
%        [  0  ]     [dtrue]
%
% where Rtrue is k-by-k upper triangular with positive diagonal.  Then we
% apply a product of random Householder reflectors on the left to obtain
% a more generic-looking (A, b) with the same geometry.
%
% INPUT:
%   n     : number of rows  (n >= k)
%   k     : number of columns
%   kind  : type of problem
%           1---> square, nonsingular system (n is reset to k)
%           2---> tall, consistent system (dtrue = 0)
%           3---> tall, inconsistent system (dtrue ~= 0)
%   dscale: scale factor for dtrue (e.g. 0.0, 1e-3, 1e-1, 1.0)
%   nref  : number of random Householder reflectors to apply on the left
%
% OUTPUT:u
%   A     : n-by-k matrix
%   b     : n-by-1 right-hand side
%   xtrue : "true" solution in the orthogonal coordinates (Rtrue * xtrue = ctrue)
%   Rtrue : k-by-k upper triangular matrix used in the construction
%   ctrue : k-by-1 vector (top part of Q'*b in the orthogonal system)
%   dtrue : (n-k)-by-1 vector (bottom part; 0 for consistent problems)
% -----------------------------------------------------------------------------------------------------

% ------------------------ Handle missing input arguments with simple defaults ----------------------
    
    if nargin < 5
        nref = 6;      % number of random reflectors
    end
    
    if nargin < 4
        dscale = 0.0;
    end
    
    if nargin < 3
        kind = 2;      % default: tall, consistent
    end

    if kind == 1
        n = k;         % Square system: force n = k
    else
        
        if n < k
            error('gen_structured_lls: require n >= k for tall systems.');
        end
    end

%====================================================================================================
    % Step 1: Build Rtrue (upper triangular with positive diagonal)
%====================================================================================================

% Start with a random k-by-k matrix
    
    Rtrue = randn(k, k);

% Zero out the entries below the diagonal so Rtrue is upper triangular
    
    for i = 1:k
        for j = 1:i-1
            Rtrue(i, j) = 0;
        end
    end

% Make the diagonal entries positive and nonzero
    
    for i = 1:k
        d = Rtrue(i, i);    % current diagonal entry
        if d < 0
            d = -d;         % make it positive
        end
        if d == 0
            d = 1;          % avoid exact zero so Rtrue is clearly nonsingular
        end
        Rtrue(i, i) = d;
    end

%====================================================================================================
    % Step 2: Choose ctrue and dtrue in the "orthogonal" basis
%====================================================================================================

% ctrue can be any vector in R^k
    
    ctrue = randn(k, 1);

% dtrue controls consistency / inconsistency of the LS problem
    
    if kind == 1    % Square system: no d-part
        
        dtrue = [];
        
    elseif kind == 2 % Tall, consistent: dtrue = 0

        dtrue = zeros(n - k, 1);
    else
       
        if n > k % Tall, inconsistent: nonzero dtrue scaled by dscale
            
            dtrue = dscale * randn(n - k, 1);
            
        else
            
            dtrue = [];
        end
    end

%     Build the simple (A0, b0) in this orthogonal coordinate system:
%             A0 = [Rtrue; 0],    b0 = [ctrue; dtrue].
    
    if n > k
        
        A0 = [Rtrue; zeros(n - k, k)];
        
    else
        
        A0 = Rtrue;
        
    end
    
    b0 = [ctrue; dtrue];

%     xtrue is the exact solution of Rtrue * x = ctrue in this system.
%         So we call our own back-substitution routine from Program 3.
    
    
    tol = 1e-12;
    xtrue = back_sub(Rtrue, ctrue, tol);

%====================================================================================================
    % Step 3: Apply random Householder reflectors on the LEFT
%====================================================================================================

    % We apply H * A0 and H * b0, where H is a product of random reflectors:
    %
    %     H = H_{nref} * ... * H_2 * H_1,
    %
    % with each H_i orthogonal of the form
    %     H_i = I - 2 u_i u_i'   (u_i is a unit vector).
    %
    % This preserves the LS geometry but produces a "generic" A and b.

    A = A0;
    b = b0;

    for t = 1:nref

       
        u = randn(n, 1);     % Random vector u in R^n

% --------------------Normalize u to have unit length: u = u / ||u||_2-------------------------------

        s = 0;
        for i = 1:n
            s = s + u(i)^2;
        end
        s = sqrt(s);

        u = u / s;   % normalize u 

        % Apply H = I - 2 u u' to A and b:
        
        %   H*A = A - 2*u*(u'*A)
        
        %   H*b = b - 2*u*(u'*b)

% ---------------------------- Compute u'*A (1-by-k row vector) -------------------------------------
        
        uTA = zeros(1, k);

        for j = 1:k              
            s_col = 0;
            for i = 1:n                          % compute the dot prod:
                s_col = s_col + u(i) * A(i, j);  % u dot A_j where A_j is the j-th col of A
            end                                  % this is the j-th comp of u'*A.
            uTA(1, j) = s_col;
        end

%--------------------------------Update A <--- A - 2*u*(u'*A)----------------------------------------
        
        for i = 1:n
            for j = 1:k
                A(i, j) = A(i, j) - 2 * u(i) * uTA(1, j);
            end
        end

%---------------------------------Compute u'*b (scalar)----------------------------------------------
        
        s_b = 0;
        for i = 1:n
            s_b = s_b + u(i) * b(i);
        end

% ------------------------------Update b <--- b - 2*u*(u'*b) ----------------------------------------
        
        for i = 1:n
            b(i) = b(i) - 2 * u(i) * s_b;
        end
    end
end


