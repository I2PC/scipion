%======================================================================== %
%
%======================================================================== %
% DEFINITION
% This function
%======================================================================== %
% INPUTS
% X
%======================================================================== %
% OUTPUTS
% Y
%======================================================================== %
% COPYRIGHT
% Authors:   laurene.donati@epfl.ch, masih.nilchian@epfl.ch
% Date:      June 2018
% Copyright: Biomedical Imaging Group, EPFL
% Citation:  Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.
%======================================================================== %

function [xVal,yVal, xGridd, yGridd]  = set2Dgrid(size, tau, scale)

% Computer center
center   = (size(1)-1)/2;

% Set X/Y vector based on refinement/scaling (if needed)
xVal = scale*(-center:tau:center);
yVal = scale*(-center:tau:center);

% Compute meshgrid
[xGridd,yGridd] = meshgrid(xVal,yVal);

end