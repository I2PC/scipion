%======================================================================== %
% [noisyIm] = addGaussianNoiseBasedOnSNR(im, SNR)
%======================================================================== %
% DEFINITION
% This function generates Gaussian noise on an image based on a desired
% output SNR. 
%======================================================================== %
% INPUTS 
% im           Noiseless image. 
% SNR          Desired output SNR level. 
%======================================================================== %
% OUTPUTS 
% noisyIm      Noisy image.  
%======================================================================== %

%======================================================================== %
% COPYRIGHT 
% Authors:   laurene.donati@epfl.ch, masih.nilchian@epfl.ch
% Date:      June 2018
% Copyright: Biomedical Imaging Group, EPFL
% Citation:  Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.
%======================================================================== %

function [noisyIm] = addGaussianNoiseBasedOnSNR(im, SNR)

    noise       = 10^(-SNR/20)*norm(im(:))/sqrt(numel(im))*randn(size(im)); 
    noisyIm     = im + noise;

end 