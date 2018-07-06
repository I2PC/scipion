%======================================================================== %
% p = KaiserBesselProjection(param, tau)
%======================================================================== %
% DEFINITION
% This function computes the 1D X-ray transform of the KBWF with a given 
% precision. (Since KBWF is isotropic, this 1D information is sufficient.)
% The close-form formula of the projection of KBWF is given by (A.7) in
% [Multidimensional digital image representation ..., Lewitt R. M., 1990]
%======================================================================== %
% INPUTS 
% param.m          Smoothness parameter of KBWF
% param.alpha      Window taper parameter of KBWF
% param.a          Support of KBWF
% tau              (positive) Fineness of the look-up table content
%======================================================================== %
% OUTPUTS 
% p          1D function. X-ray transform of the KBWF along one radius. 
%======================================================================== %

%======================================================================== %
% COPYRIGHT 
% Authors:   masih.nilchian@epfl.ch, laurene.donati@epfl.ch 
% Date:      June 2018
% Copyright: Biomedical Imaging Group, EPFL
% Citation:  Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.
%======================================================================== %

function p = KaiserBesselProjection1D(param, tau)

% Get KBWF information 
m      = param.m; 
alpha  = param.alpha; 
a      = param.a; 

% Set fineness of look-up table  
s = 0:tau:a; 

% Compute 1D Projection of KBWF 
tmp =  sqrt(1 -(s/a).^2);
p   =  a ./ besseli(m, alpha) .* sqrt(2*pi/alpha) .* tmp.^(m+0.5) .* besseli(m+0.5, alpha*tmp);

end 

