%======================================================================== %
% ctf =  create2Dctf(sizeCTF, paramCTF)
%======================================================================== %
% DEFINITION
% This function generates a 2D contrast transfer function (CTF).
%======================================================================== %
% INPUTS 
% sizeCTF            Size of 2D CTF (e.g., size of projections). 
% paramCTF.sigma     Standard deviation of CTF image. 
% paramCTF.type      Type of CTF shape ('gaussian' or 'sinc') 
%======================================================================== %
% OUTPUTS 
% ctf               2D CTF image.  
%======================================================================== %

%======================================================================== %
% COPYRIGHT 
% Authors:   laurene.donati@epfl.ch, masih.nilchian@epfl.ch
% Date:      June 2018
% Copyright: Biomedical Imaging Group, EPFL
% Citation:  Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.
%======================================================================== %


function ctf =  create2Dctf(sizeCTF, paramCTF)

% Get ctf parameters
sigma = paramCTF.sigma; 
type  = paramCTF.type;

% Set grid information
supp   = floor((1+sizeCTF(1))/2);
[X,Y]  = meshgrid(linspace(-supp,supp,1));

% Sample continuous CTF on discrete grid 
if strcmp(type, 'gaussian')
    ctf    = exp(-(X.^2+Y.^2)/(2*sigma.^2));
elseif strcmp(type, 'sinc')
    ctf    = sinc(-(X.^2+Y.^2)/(2*sigma.^2));
end

end