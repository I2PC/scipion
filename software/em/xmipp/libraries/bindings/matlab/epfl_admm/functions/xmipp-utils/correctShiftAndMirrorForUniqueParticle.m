% laurene.donati@epfl.ch


% laurene.donati@epfl.ch
% This function corrects the shift of the projections through spline interpolation
% The shift information are output from Scipion.
% NOTE: Currently works only for square projection images.

function proj = correctShiftAndMirrorForUniqueParticle(proj,shiftX,shiftY,flipIndex)


% Padd projections to handle shift at border
maxShift    = ceil(max(abs([shiftX; shiftY])));
paddedProj = padarray(proj,[maxShift,maxShift],0,'both');

% Get size info
Nx  = size(proj,1);
Nxp = size(paddedProj,1);

% Get centered projection meshgrid
cP        = (Nxp-1)/2;
xValP     = -cP:cP;
yValP     = -cP:cP;
[pGx,pGy] = meshgrid(xValP,yValP);

% Get interpolation center
cI = (Nx-1)/2;

% Handle shifts and projection mirrors (TO TEST!!!)
if ~flipIndex
    xValI     = (-cI:cI) - shiftX;
    yValI     = (-cI:cI) - shiftY;
else
    xValI     = -(-cI:cI) + shiftX;
    yValI     = (-cI:cI) - shiftY;
end

% Create interpolation meshgrid
[iGx,iGy] = meshgrid(xValI,yValI);

% Interpolate with sinc
proj = interp2(pGx,pGy,paddedProj,iGx,iGy, 'sinc');


end

