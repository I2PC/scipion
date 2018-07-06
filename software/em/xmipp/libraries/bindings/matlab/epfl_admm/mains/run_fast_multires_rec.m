%=========================================================================
% run_fast_multires_rec.m
%=========================================================================
% DEFINITION 
% This script runs our fast multiresolution reconstruction algorithm for
% cryo-EM on synthetic data. Projections are computed and a volume is
% is reconstructed at a given scale. 
%=========================================================================
% MAIN INPUTS 
% -------------- for simulation 
% volName      (String) Name of 3D ground-truth volume.
% sizeIndex    (Scalar) Size along 1D of 3D ground-truth volume.
% nberProjs    (Scalar) Number of simulated projection directions. 
% snrProjs     (Scalar) SNR of simulated projection images. 
% ctf.use      (Boolean) Simulates CTF effect on projection images. 
% ctf.type     (String) Type of simulated CTF ('gaussian', 'sinc')
% ctf.sigma    (Double) Standard deviation of simulated CTF
% tauProj      (Double) 
% -------------- for reconstruction 
% scale        (Scalar) Scale of the reconstructed volume.
% tauHtb       (Double) 
% tauKernel    (Double) 
% alphaRec     (Scalar) Regularization parameter. 
% nbItADMM     (Scalar) Number of outer ADMM iterations. 
% nbItCG       (Scalar) Number of inner ADMM iterations.  
%=========================================================================
% MAIN OUTPUTS 
% rec          (Volume) Reconstructed 3D volume. 
%=========================================================================
% NOTES 
% - Currently valid only for cubic volumes. 
%=========================================================================

%=========================================================================
% Authors: laurene.donati@epfl.ch, masih.nilchian@epfl.ch
% Date: June 2018
% Copyright: Biomedical Imaging Group, EPFL
% Citation: Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.
%=========================================================================

clear all; clc; close all;

% Set default paths 
indComp = 4;
switch indComp
    case 1 % Personal computer (Mac)
        mainPath = '/Users/laurenedonati/Dropbox/multires-cryo';
    case 2 % Personal laptop (Mac)
        mainPath = '/Users/laurene/Dropbox/multires-cryo';
    case 3 % Masih's laptop
        mainPath = '/Users/masih.nilchian/Dropbox/multires-cryo';     
    case 4 % New Linux laptop 
        mainPath = '/home/big/Desktop/multires-cryo';  
end 
addpath(genpath(mainPath)); cd([mainPath '/functions']);


%% ========================================================================
%                  PART 0.  PARAMETERS
% =========================================================================

% Simulation parameters
volName       = 'betagal';
sizeIndex     = 64;
nberProjs     = 20;
snrProjs      = 100; 

% CTF parameters 
ctf.use = false;
if ctf.use
    ctf.type   = 'gaussian';
    ctf.sigma  = 5;
end

% Reconstruction parameters
scale         = 2;
alphaRec      = 1e-7;
nbItADMM      = 30;
nbItCG        = 5;

% Refinement paramaters for LT
tauProj     = 0.01;
tauHtb      = scale*0.25; 
tauKernel   = 0.01;

% Booleans
DISPLAY_FIGS  = true;


%% ========================================================================
%                  PART 1.   INITIALISATION
% =========================================================================

% KBWF parameters (used for projection*)
kbwfProj    = struct('a', 2, 'alpha', 10.8, 'm',2);

% KBWF parameters (used for reconstruction*)
kbwfRec     = struct('a', 4, 'alpha', 19, 'm', 2);

% Set gt size
sizeVol     = [sizeIndex; sizeIndex; sizeIndex];

% Diverse paths
gtPath      = [mainPath '/data/volumes/' volName '/' volName ...
                    '-' num2str(sizeIndex) '.mat'];
KBWFcoeffPath  = setKBWFCoeffPath(mainPath,volName,sizeIndex);


%% ========================================================================
%                  PART 2.   SIMULATION
% =========================================================================

% Load GT (for metric computation) and KBWF coeffs
vol       = struct2array(load(gtPath));
ckVol     = struct2array(load(KBWFcoeffPath));

% Set equi-distributed projection directions
[rot,tilt,psi] = generateEquidistributedAngles(nberProjs);
[r1,r2]        = setProjectionPlanes(rot,tilt,psi, false);

% Compute look-up table with KBWF projection
ltProj1D = KaiserBesselProjection1D(kbwfProj, tauProj);

% Compute projection measurements
disp('-- Compute projections.');
projs    = projection3DMT(ckVol, r1, r2, [1,1,1], [1,1], ltProj1D, kbwfProj.a);

% Convolve projections with CTF
if ctf.use == true
    ctf.image = create2Dctf(size(projs), ctf);
    projs     = convolveProjsWithCTF(projs,ctf);
end

% Add Gaussian noise on projections
[projs]  = addGaussianNoiseBasedOnSNR(projs, snrProjs); 


%% ========================================================================
%                  PART 3.   RECONSTRUCTION
% =========================================================================
 
rec = reconstruct_multires_ADMM(projs, r1, r2, scale, alphaRec, nbItADMM, nbItCG, sizeIndex);

% % Set reconstruction size (based on scale)
% sizeRec    = floor((sizeVol-1)/scale+1);
% 
% % Combut Htb (fast version)
% Htb        = computeHtbfast(projs, kbwfRec, tauHtb, r1, r2, scale, sizeRec, 0);
% 
% % Compute kernel and define HtH (with kernel)
% kernel     = computeKernelHtH3D(kbwfRec,tauKernel, r1, r2, scale, sizeRec, ctf);
% HTH        = @(X) applyKernel3D(kernel,X);
% 
% % Set regularization (prox operator)
% prox       = @(X,alphaRec) softThresholdingI3D(X,alphaRec);
% 
% % Set projection constraints
% projConstr = @(X) positiveProj(X);
% 
% % Run ADMM reconstruction 
% tic
% [ckRec,~]  = applyADMM(Htb,'ATA',HTH,'verbose', true,'drawflag',false, ...
%         'alpha',alphaRec,'lambda',0,'tol',1e-20,'maxOutIter', nbItADMM, ...
%     'maxInIter',nbItCG, 'mu',1e1*alphaRec,'doProjection',projConstr,'prox',prox);
% 
% % Re-expand reconstructed ck in fine grid
% rec        = expandScaledKBWFcoeffstoFineGrid_3D(ckRec, scale, kbwfRec, sizeVol);
% 
% % Display reconstruction time 
% disp(['Time to run ADMM and re-expand: ' num2str(toc) ' sec'])


%% ========================================================================
%                  PART 4.   PROCESS RESULTS
% =========================================================================

% Display orthoslice, FSC curves and profile lines 
if DISPLAY_FIGS
   display_result_multires_rec(vol,rec,sizeIndex);
end


% =========================================================================
% FOOTNOTES 
% *: To reduce the risk of inverse crime. 

