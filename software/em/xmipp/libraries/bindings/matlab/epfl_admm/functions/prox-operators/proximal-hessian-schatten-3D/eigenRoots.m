%Function eigenRoots finds the eigenvalues of a 3x3 symmetric matrix.
%
%Matlab Usage: E=eigenRoots(Y);
%
% Function eigenRoots can handle multiple matrix inputs.
% Y: Y is a multidimensional matrix where its first dimension should be of
% size 6, corresponding to the 6 unique elements of the 3x3 symmetric
% matrix, i.e., 
% Y(:,k) = [Y(1,1,k) Y(1,2,k) Y(1,3,k) Y(2,2,k) Y(2,3,k) Y(3,3,k)]
% 