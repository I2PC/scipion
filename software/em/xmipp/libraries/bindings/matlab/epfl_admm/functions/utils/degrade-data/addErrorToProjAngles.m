
% laurene.donati@epfl.ch

function [noisyAngles] = addErrorToProjAngles(angles, var)

[dim1, NberAngles] = size(angles);
errorArray = sqrt(var)*randn(dim1, NberAngles);
noisyAngles = double(angles) + errorArray;

end 
 