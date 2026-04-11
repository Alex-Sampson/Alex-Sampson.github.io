% Below are the test problems for Task 1
% The given inputs are all in-order coordinate form of sparse matrices
% The vectors to multiply in each case are given in dense non-sparse  form

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrix 1
AA = [10,-1,3,11,5,12,7,13,9,2,3];
JR = [1,1,2,2,3,3,4,4,5,5,5];
JC = [1,5,1,2,2,3,3,4,4,1,5]; %Column out of order at the end

%Vector 1 to Multiply Against Matrix 1
v = [1; 0; -1; 0; 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrix 2
AA2 = [1, 5, -2, 4, 1, 2, 1, 2, 3];
JR2 = [1, 1, 2, 2, 3, 3, 4, 4, 4];
JC2 = [1, 4, 1, 2, 1, 4, 2, 3, 4];

%Vector 2 to Multiply Against Matrix 2
v2 = [3; 5; -1; -1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrix 3
AA3 = [1, 2, 3, 2, -3, 1, 2, 3, -1, 1, 2, 1, -1, 4, 1, 1, -4, 5];
JR3 = [1, 1, 1, 2, 2, 3, 3, 4, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8];
JC3 = [1, 3, 5, 2, 6, 3, 5, 4, 2, 5, 3, 4, 5, 6, 2, 7, 4, 8];

%Vector 3 to Multiply Against Matrix 3
v3 = [-2; 2; -2; 2; -2; 2; -2; 2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
