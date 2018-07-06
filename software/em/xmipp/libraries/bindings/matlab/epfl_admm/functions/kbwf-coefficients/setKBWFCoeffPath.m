%======================================================================== %
% coeffPath = setKBWFCoeffPath(mainPath, volName,volSize)
%======================================================================== %
% DEFINITION
% This function returns the path to the KBWF coefficients of synthetic 
% volumes based on their sizes. 
%======================================================================== %
% INPUTS 
% mainPath     Path leading to main cryo-multires folder 
% volName      Name of volume (e.g., 'betagal')
% volSize      Size index of cubic volume (e.g., 128)
%======================================================================== %
% OUTPUTS 
% coeffPath  Path to the KBWF coefficients
%======================================================================== %
% COPYRIGHT 
% Authors:   laurene.donati@epfl.ch, masih.nilchian@epfl.ch
% Date:      June 2018
% Copyright: Biomedical Imaging Group, EPFL
% Citation:  Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.
%======================================================================== %


function coeffPath = setKBWFCoeffPath(mainPath, volName, volSize)

% Set number of ADMM iterations for computing KBWF coefficients
if volSize == 256
    nberKBWFiterADMM = 100;
elseif volSize == 128
    nberKBWFiterADMM = 500;
elseif volSize == 64
    nberKBWFiterADMM = 1000;
end      
       
% Path to KBWF coefficients 
coeffPath    = [mainPath '/data/kbwf-coeffs/' volName '/' volName ...
                '-' num2str(volSize) '-kbwfCoeffs-It' num2str(nberKBWFiterADMM) '.mat'];

end