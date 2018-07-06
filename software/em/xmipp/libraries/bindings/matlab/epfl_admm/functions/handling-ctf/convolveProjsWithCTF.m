%======================================================================== %
% projs   = convolveProjsWithCTF(projs, paramCTF)
%======================================================================== %
% DEFINITION
% This function convolves a 2D CTF with a set of projection images.
%======================================================================== %
% INPUTS 
% projs             Set of 2D projection measurements. 
% paramCTF.use      (Boolean) Set whether CTF effect is considered. 
% paramCTF.image    2D CTF image. 
%======================================================================== %
% OUTPUTS 
% ctf               Set of convolved 2D projection measurements.  
%======================================================================== %
% NOTES 
% - Considers a unique CTF for all projs as of now. 
%========================================================================

%======================================================================== %
% COPYRIGHT 
% Authors:   laurene.donati@epfl.ch, masih.nilchian@epfl.ch
% Date:      June 2018
% Copyright: Biomedical Imaging Group, EPFL
% Citation:  Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.
%======================================================================== %

function projs         = convolveProjsWithCTF(projs,paramCTF)

if paramCTF.use == true
    image = paramCTF.image; 
    parfor i = 1: size(projs,3)
        projs(:,:,i)  =  real(ifft2(fft2(squeeze(projs(:,:,i))).*fftshift(image)));
    end
end

end