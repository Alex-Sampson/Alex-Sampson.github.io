%
% This is a test driver for the
% testing the multiplication subroutines
% (and their different versions)
% 
% Lv-->z,  Uv-->z, LU-->A
%
%
% The actions of this routine are controlled
% by the parameters set at the beginning.
%
% Tests are run for a specifed subroutine
% for a range of values of the dimension n
% (the outer loop of code)
%
% For each value of n matrices L and U
% are generated along with a vector v.
% (the main inner loop of code)
% This can be done based on a custom 
% specification of these matrices as a function
% of n (which the user must place in the code)
% or a random number generator based matrix
% with some structural constraints imposed.
% Depending on the subroutine tested,
% L and v, or U and v, or L and U will be 
% used first as Matlab objects to use
% Matlab's library and commands to determine
% the "true" answer for comparison to the
% computed answer.
% The matrices and vectors required are then
% placed in the data structures as required
% by the assignment, the test is run,
% the error information collected.
%
% During this loop, if desired, detailed
% information on the results of each test problem
% can be output. Or, a summary (mean and max) of the
% relative error can be output for each set of tests
% with the same n. A third output option collects the
% information and plots the means and maximums of each set of
% tests with the same n vs the value of n.
% In all cases at the end the maximum value of the relative
% error over all tests is printed.
%

%
% The arrays M and Acomp,  are not n by n in order
% to have a code similar to that used in an
% when a range of problems would use the same
% data structure or when other information
% would be stored in M and Acomp in addition
% to L, U, and L*U. The same is true for the
% vectors v, zcomp, and ztrue.
%


%
% These set the dimensions of the  data structures.
% All problems tested should have n smaller than the smallest of these.
Mrdim=nmax;
Mcdim=nmax;
Acomprdim=nmax;
Acompcdim=nmax;
zdim=nmax;
vdim=nmax;

%
% set range of matrix dimensions and number
% of tests per dimension
%
nmin=30;
nmax=900;
ndelta=10;
ntests=5;

%
% Indicate the subroutine to be tested by commenting out those not desired.
%
% routinenum = 0; % L*U row-oriented outer product algorithm          <<-- Use
% routinenum = 1; % L*U row-oriented inner product algorithm          <<-- Wont Use
% routinenum = 2; % L*U column-oriented vector triad  algorithm       <<-- Wont Use
% routinenum = 3; % L*v row-oriented inner product algorithm          <<-- Use
% routinenum = 4; % L*v column-oriented vector triad  algorithm       <<-- Use
% routinenum = 5; % U*v row-oriented inner product algorithm          <<-- Given
% routinenum = 6; % U*v column-oriented vector triad  algorithm       <<-- Given
% routinenum = 7; % L*v compressed-column (1D) vector triad algorithm <<-- Added
% routinenum = 8; % L*v banded multiplication algorithm               <<-- Added
%
% Set the manner of matrix and vector generation by commenting out those not desired.
%
%problemtype = 1; % User custom manner. Add the code at the appropirate place below.
%problemtype = 2; % Uniform random numbers on the interval [dmin,dmax].
                 % Note that this interval's endpoints can be set to include both signs
                 % or have a fixed one based on the signs of the endpoints.
                 % this does not affect the signs of the 1's on the diagonal of L.
%dmin=-10.0;
%dmax=10.0;
problemtype = 3; % Normal random numbers with specified mean and standard deviation.
tempmean=0.0;
tempstdev=500.0;

%-----------------------------------------------------------------------------------------
% Set output type by commenting out those not desired.
%
 %outputtype = 1; % norm information listed for every test problem
 %outputtype = 2; % only mean/max norm information listed for each set with fixed n
 outputtype = 3; % mean/max norm vs n plotted (the arrays are available for 
                % examination after driver runs in maxforthisn and meanforthisn arrays
%-----------------------------------------------------------------------------------------

%
% There are no parameters for the user to set below this point
% except the problemtype = 1 code to set the the manner in which
% the problem is generated, i.e., the set of numbers from which
% L and U are selected and the vector v.
%

%
% initialize max relative error for all  tests of all sizes
%
maxrelerror = 0.0;

%
% Output Label of Tests to the Matlab Terminal
%
if routinenum == 0
  fprintf(1,'\n**************************************** \n')
  fprintf(1,'A = L*U Outer Product Version Tests for sizes n = %g to %g by %g \n',nmin,nmax,ndelta)
  fprintf(1,'**************************************** \n')
elseif routinenum == 1
  fprintf(1,'\n**************************************** \n')
  fprintf(1,'A = L*U Row Version Tests for sizes n = %g to %g by %g \n',nmin,nmax,ndelta)
  fprintf(1,'**************************************** \n')
elseif routinenum == 2 
  fprintf(1,'\n**************************************** \n')
  fprintf(1,'A = L*U Column Version Tests for sizes n = %g to %g by %g \n',nmin,nmax,ndelta)
  fprintf(1,'**************************************** \n')
elseif routinenum == 3
  fprintf(1,'\n**************************************** \n')
  fprintf(1,'z = L*v Row  Version Tests for sizes n = %g to %g by %g \n',nmin,nmax,ndelta)
  fprintf(1,'**************************************** \n')
elseif routinenum == 4
  fprintf(1,'\n**************************************** \n')
  fprintf(1,'z = L*v Column  Version Tests for sizes n = %g to %g by %g \n',nmin,nmax,ndelta)
  fprintf(1,'**************************************** \n')
elseif routinenum == 5
  fprintf(1,'\n**************************************** \n')
  fprintf(1,'z = U*v Row  Version Tests for sizes n = %g to %g by %g \n',nmin,nmax,ndelta)
  fprintf(1,'**************************************** \n')
elseif routinenum == 6
  fprintf(1,'\n**************************************** \n')
  fprintf(1,'z = U*v Column  Version Tests for sizes n = %g to %g by %g \n',nmin,nmax,ndelta)
  fprintf(1,'**************************************** \n')
elseif routinenum == 7
  fprintf(1,'\n**************************************** \n')
  fprintf(1,'z = L*v Compressed Column Version Tests for sizes n = %g to %g by %g \n',nmin,nmax,ndelta)
  fprintf(1,'**************************************** \n')
elseif routinenum == 8
  fprintf(1,'\n**************************************** \n')
  fprintf(1,'z = L*v Banded Version Tests for sizes n = %g to %g by %g \n',nmin,nmax,ndelta)
  fprintf(1,'**************************************** \n')
else
  fprintf(1,'\n**************************************** \n')
  fprintf(1,'Undefined version number   = %g \n',routinenum)
  fprintf(1,'**************************************** \n')
  return; % terminate the tester and return to command prompt
end

%
% Initialize error saving vectors  
% 
%
Error2Norm = zeros(ntests,1);
Atrue2Norm = zeros(ntests,1);
RelError =  zeros(ntests,1);
ztrue2Norm = zeros(ntests,1);

% Initialize arrays to collect mean/max norm for every set with fixed n
nexperiments=size([nmin:ndelta:nmax]',1);
meanforthisn = zeros(nexperiments,1);
maxforthisn =  zeros(nexperiments,1);
experimentindex = 0;

%Test loop over a range of dimensions of L, U, and v
for n = nmin:ndelta:nmax

  experimentindex = experimentindex + 1;
  if (outputtype == 1) || (outputtype == 2)
    fprintf(1,'\n********************************** \n')
    fprintf(1,'Tests for size n = %g \n',n)
    fprintf(1,'********************************** \n')
  end

% loop over fixed n for ntests random problems
  for itest=1:ntests
% reinitialize data structures for a new problem with size n
    M=zeros(Mrdim,Mcdim);
    Acomp=zeros(Acomprdim,Acompcdim);
    zcomp=zeros(zdim,1);
    ztrue=zeros(zdim,1);
    v=zeros(vdim,1);
%
% Generate full expanded L and U
% along with a vector v
% and place them each of them in their
% associated work array.
%
% Note that since this is testing mulitplication routines
% there is no need to guarantee that the diagonal elements
% of U are nonzero.
%

% Save current state of random number generator.
% This can be used for reproducing the same sequence of
% test problems if needed. That is not done in this
% form of the code but it is often useful for debugging.
% the rng state can be restored with the statement rng(statesave);
% See the Matlab documentation for more details.
  statesave=rng;

% Select the type of random numbers 
% and set the title of the error plot for outputtype = 3
% by commenting/uncommenting the appropriate lines

  if problemtype == 1
%  
    outputtitle = "PUT YOUR CUSTOM TITLE HERE";
    Temp = ones(n,n); % . Replace with yours.
    v(1:n)=randn(n,1); % This is where you would set v all 1s
  elseif problemtype == 2
% standard uniform on (0,1) shifted to (dmin,dmax) set earlier
    strdmin=string(dmin);
    strdmax=string(dmax);
    temptitle = ["uniform distribution with shifted to (" strdmin "," strdmax ")"];
    outputtitle = join(temptitle);
    Temp = dmin.*ones(n,n) + (dmax-dmin).*(rand(n,n)); 
    v(1:n) = dmin.*ones(n,1) + (dmax-dmin).*(rand(n,1)); 
  elseif problemtype == 3
%  normal distribution with  mean and standard deviation set
%  to tempmean and tempstdev set earlier
    Temp  = tempstdev.*randn(n,n) + tempmean.*ones(n,n);
    v(1:n) = tempstdev.*randn(n,1) + tempmean.*ones(n,1);
    strmean=string(tempmean);
    strstdev=string(tempstdev);
    temptitle = ["normal element distribution with mean/deviation " strmean "and" strstdev];
    outputtitle = join(temptitle);
  else
    fprintf(1,'\n**************************************** \n')
    fprintf(1,'Undefined problem number   = %g \n',routinenum)
    fprintf(1,'**************************************** \n')
    return; % terminate the tester and return to command prompt
  end

  L = eye(n,n)+ tril(Temp,-1); % L is strict lower part of 
                               % the array of random numbers
                               % and a diagonal of all 1
  U = triu(Temp);  % U is strict lower part of 
                   % the array of random numbers


%
%  Generate true product of L and U using the built-in
%  Matlab matrix product operator and store the result
%  a work array Atrue not involved in your subroutines.
%

  if (routinenum == 0) || (routinenum == 1) || (routinenum == 2) 
    Atrue = L*U;
  elseif (routinenum == 3) || (routinenum == 4)|| (routinenum == 7)
    ztrue(1:n) = L*v(1:n);
    
  elseif routinenum == 8

%----Code for extracting bands Lb1 and Lb2 from our 'banded' unit lower triangular L:----
Lb1 = diag(L, -1);       % length n-1
Lb2 = diag(L, -2);       % length n-2
LB=eye(n);               % initializes our banded matrix
for i=1:size(Lb1,1)
    LB(i+1,i)=Lb1(i);    % populates first band
end
for j=1:size(Lb2,1)
    LB(j+2,j)=Lb2(j);    % populates second band
end
%-----------------------------------------------------------------------------------------
    ztrue(1:n) = LB * v(1:n); % true calculation for Banded Matrix vector product.

  elseif (routinenum == 5) || (routinenum == 6)
    ztrue(1:n) = U*v(1:n);
  else
    fprintf(1,'\n**************************************** \n')
    fprintf(1,'Undefined routine number   = %g \n',routinenum)
    fprintf(1,'**************************************** \n')
    return; % terminate the tester and return to command prompt
  end

%
% Put L and U into the  simple compressed
% storage of L and U % in the array M.
% Note the 0s of L and U are not stored and
% the 1's of the diagonal of L are note stored in the
% work array M(mrdim,mcdim).
% The vectors v, zcomp, and ztrue are already in work arrays.
%
%     First a dense write of U, i.e., this includes elements of
%     the U array that are not actually in the mathematical matrix U
%     then overwrite 0s in U and M  with lower part of L.
%

%---------Code for constructing M from the random matrices L and U------------------

  M(1:n,1:n) = U;
  for i = 1:n-1
     M(i+1:n,i)= L(i+1:n,i);
  end
%-----------------------------------------------------------------------------------

%---------------Code for writing L as a compressed column array:--------------------
  Lc = zeros(n*(n-1)/2, 1); % initialize an array of zeros of appropriate length 
  k = 1; % Indexing counter
  for j = 1:n-1           % columns from the first to the 2nd to last
    for i = j+1:n         % rows below the diagonal in column j
        Lc(k) = M(i,j);   % Assigns each entry of Lc one column of M at a time
        k = k + 1;        % Increases our index step after each entry assignment.
    end
  end
%-----------------------------------------------------------------------------------


% Run test number itest for the current value of n
%
  if (routinenum == 0) 
    [Acomp]=LUmult_outer(M, Acomp, n);
    ErrorMatrix = Acomp(1:n,1:n) - Atrue(1:n,1:n);
%
%   Compute the norms for this problem and save it with the others with this n
%
    Error2Norm(itest) = norm(ErrorMatrix,2);
    Atrue2Norm(itest) = norm(Atrue,2);
    RelError(itest) = Error2Norm(itest)/Atrue2Norm(itest);
%
%   update running maximum relative error over test problems for all n
%
    maxrelerror = max(maxrelerror,RelError(itest));
  elseif (routinenum == 1) 
    [Acomp]=LUmult_row(M, Acomp, n);
    ErrorMatrix = Acomp(1:n,1:n) - Atrue(1:n,1:n);
%
%   Compute the norms for this problem and save it with the others with this n
%
    Error2Norm(itest) = norm(ErrorMatrix,2);
    Atrue2Norm(itest) = norm(Atrue,2);
    RelError(itest) = Error2Norm(itest)/Atrue2Norm(itest);
%
%   update running maximum relative error over test problems for all n
%
    maxrelerror = max(maxrelerror,RelError(itest));
  elseif (routinenum == 2) 
    [Acomp]=LUmult_col(M, Acomp, n);
    ErrorMatrix = Acomp(1:n,1:n) - Atrue(1:n,1:n);
%
%   Compute the norms for this problem and save  with the others with this n
%
    Error2Norm(itest) = norm(ErrorMatrix,2);
    Atrue2Norm(itest) = norm(Atrue,2);
    RelError(itest) = Error2Norm(itest)/Atrue2Norm(itest);
%
%   update running maximum relative error over test problems for all n
    maxrelerror = max(maxrelerror,RelError(itest));
  elseif (routinenum == 3)
    [zcomp]=Lvmult_row(M, v, zcomp,n);
    ErrorVector = zcomp(1:n) - ztrue(1:n);
%
%   Compute the norms for this problem and save with the others with this n
%
    Error2Norm(itest) = norm(ErrorVector,2);
    ztrue2Norm(itest) = norm(ztrue,2);  % ztrue is an n element vector
    RelError(itest) = Error2Norm(itest)/ztrue2Norm(itest);

% update running maximum relative error over test problems for all n
    maxrelerror = max(maxrelerror,RelError(itest));
  elseif (routinenum == 4)
    [zcomp]=Lvmult_col(M, v, zcomp,n);
    ErrorVector = zcomp(1:n) - ztrue(1:n);
%
%   Compute the norms for this problem and save with the others with this n
%
    Error2Norm(itest) = norm(ErrorVector,2);
    ztrue2Norm(itest) = norm(ztrue,2);  % ztrue is an n element vector
    RelError(itest) = Error2Norm(itest)/ztrue2Norm(itest);
% update running maximum relative error over test problems for all n
    maxrelerror = max(maxrelerror,RelError(itest));
  elseif (routinenum == 5) 
    [zcomp]=Uvmult_row(M, v, zcomp,n);
    ErrorVector = zcomp(1:n) - ztrue(1:n);
%
%   Compute the norms for this problem and save with the others with this n
%
    Error2Norm(itest) = norm(ErrorVector,2);
    ztrue2Norm(itest) = norm(ztrue,2);  % ztrue is an n element vector
    RelError(itest) = Error2Norm(itest)/ztrue2Norm(itest);
%   update running maximum relative error over test problems for all n
    maxrelerror = max(maxrelerror,RelError(itest));
  elseif (routinenum == 6)
    [zcomp]=Uvmult_col(M, v, zcomp,n);
    ErrorVector = zcomp(1:n) - ztrue(1:n);
%
%   Compute the norms for this problem and save with the others with this n
%
    Error2Norm(itest) = norm(ErrorVector,2);
    ztrue2Norm(itest) = norm(ztrue,2);  % ztrue is an n element vector
    RelError(itest) = Error2Norm(itest)/ztrue2Norm(itest);
% update running maximum relative error over test problems for all n
    maxrelerror = max(maxrelerror,RelError(itest));

  elseif (routinenum == 7)
    [zcomp]=Lvmult_compcol(Lc, v, zcomp,n);
    ErrorVector = zcomp(1:n) - ztrue(1:n);
%
%   Compute the norms for this problem and save with the others with this n
%
    Error2Norm(itest) = norm(ErrorVector,2);
    ztrue2Norm(itest) = norm(ztrue,2);  % ztrue is an n element vector
    RelError(itest) = Error2Norm(itest)/ztrue2Norm(itest);
% update running maximum relative error over test problems for all n
    maxrelerror = max(maxrelerror,RelError(itest));

    elseif (routinenum == 8)
    [zcomp]=Lvmult_banded(Lb1, Lb2, v, zcomp, n);
    ErrorVector = zcomp(1:n) - ztrue(1:n);
%
%   Compute the norms for this problem and save with the others with this n
%
    Error2Norm(itest) = norm(ErrorVector,2);
    ztrue2Norm(itest) = norm(ztrue,2);  % ztrue is an n element vector
    RelError(itest) = Error2Norm(itest)/ztrue2Norm(itest);
% update running maximum relative error over test problems for all n
    maxrelerror = max(maxrelerror,RelError(itest));
  
  else
    fprintf(1,'\n**************************************** \n')
    fprintf(1,'Undefined routine number   = %g \n',routinenum)
    fprintf(1,'**************************************** \n')
    return; % terminate the tester and return to command prompt
  end

  end % end of itest loop

if outputtype == 1
  for iouts = 1:ntests
   if (routinenum == 0) || (routinenum == 1) || (routinenum == 2)
    fprintf(1,'n = %g Error 2-norm = %5.2e  Atrue 2-norm = %5.2e Relative Error = %5.2e\n',n,Error2Norm(iouts), Atrue2Norm(iouts), RelError(iouts))
   else
    fprintf(1,'n = %g Error 2-norm = %5.2e  ztrue 2-norm = %5.2e Relative Error = %5.2e\n',n,Error2Norm(iouts), ztrue2Norm(iouts), RelError(iouts))
   end % end of if
  end % output loop
end %  end outputtype ==1 if
meanforthisn(experimentindex) = mean(RelError(1:ntests));
maxforthisn(experimentindex) =  max(RelError(1:ntests));

if outputtype == 2
 fprintf(1,'For %g problems with n = %g  Mean Relative Error = %5.2e Maximum Relative Error = %5.2e\n',ntests,n, meanforthisn(experimentindex), maxforthisn(experimentindex))
end

end % end of n loop

% print overall maximum relative error
fprintf(1,'\n\n\nMaximum Relative Error over all tests of all sizes = %5.2e\n',maxrelerror);

if outputtype == 3
% protect against -infinity in log vectors
tempmean(1:nexperiments)=max(-17,log10(meanforthisn(1:nexperiments)));
tempmax(1:nexperiments)=max(-17,log10(maxforthisn(1:nexperiments)));
plot([nmin:ndelta:nmax],tempmean(1:nexperiments),'-.k', 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 4) % plot 
hold on
plot([nmin:ndelta:nmax],tempmax(1:nexperiments),'-c','LineWidth', 1)  % plot 
xlabel('dimension of matrices in set')
ylabel('log base 10 of relative error max and mean per set')
legend('Mean RelError','Max RelError','Location','best','FontSize',11)
title(outputtitle);

% debug printing of mean and max data
% fprintf(1,'\n\n\nMean Relative Errors vs n = %5.2e\n',meanforthisn(1:nexperiments));
% fprintf(1,'\n\n\nMaximum  Relative Errors vs n = %5.2e\n',maxforthisn(1:nexperiments));
hold off
end