
% laurene.donati@epfl.ch
% TO TEST: symLR = '/home/big/LR.txt';
% NOTE: ADD AS FIRST INPUT OF Rs AN IDENTITY MATRIX TO CONSIDER
% UNSYMMETRIZED PROJECTION
% No symmetry is handled if symLR is is an empty text file 

function [Rs, nSymm]  = readSymmetryMatricesXmipp(symLR)

% Read text file
fileID = fopen(symLR,'r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
fclose(fileID);

if isempty(A)
    % If no symmetry, return identity matrix 
    nSymm = 1;
    Rs    = zeros(3,3,nSymm);
    Rs(:,:,1) = eye(3);
else
    % Sort information in matrices
    nSymm = (length(A)/(16*2))+1;
    
    % Initialise Rs
    Rs = zeros(3,3,nSymm);
    
    % Set first Rs as identity (unsymmetrized volume)
    Rs(:,:,1) = eye(3);
    
    % Read data (only R matrices)
    for i = 1:(nSymm-1)
        N = (2*i-1)*16;
        list = A((N+1):(N+16));
        temp = reshape(list, [4,4]);
        Rs(:,:,i+1) = temp(1:3,1:3);
    end
end

end