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

function [xVal, yVal, zVal, xGridd, yGridd, zGridd]  = set3Dgrid(size, tau, scale)

% Computer center
center   = (size(1)+1)/2;

% Compute step size
stepSize =  tau*scale;

% Set X/Y vector based on refinement/scaling (if needed)
xVal = -center:stepSize:center;
yVal = -center:stepSize:center;
zVal = -center:stepSize:center;

% Compute meshgrid
[xGridd, yGridd, zGridd] = meshgrid(xVal,yVal, zVal);

end