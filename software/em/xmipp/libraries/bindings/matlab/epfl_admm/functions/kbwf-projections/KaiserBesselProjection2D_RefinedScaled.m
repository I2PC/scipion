%======================================================================== %
% kbwf = KaiserBesselProjection2D_RefinedScaled(kbwf,tau, scale,sizeProj)
%======================================================================== %
% DEFINITION
% This function ...
%======================================================================== %
% INPUTS
% kbwf          
% tau 
% scale
% sizeProj 
%======================================================================== %
% OUTPUTS
% kbwf       Def 
%======================================================================== %

%======================================================================== %
% COPYRIGHT
% Authors:   laurene.donati@epfl.ch, masih.nilchian@epfl.ch
% Date:      June 2018
% Copyright: Biomedical Imaging Group, EPFL
% Citation:  Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.
%======================================================================== %


function kbwf   = KaiserBesselProjection2D_RefinedScaled(kbwf,tau, scale,sizeProj)

% Image sizes 
Nx  = sizeProj(1); 
Ny  = sizeProj(2); 

% Refined grid 
center     = (Nx-1)/2;
xValRef    = -center:tau:center; 

% Compute projection 
kbwf     = abs(scale)* KaiserBesselProjection2D(kbwf, xValRef/scale);


end 