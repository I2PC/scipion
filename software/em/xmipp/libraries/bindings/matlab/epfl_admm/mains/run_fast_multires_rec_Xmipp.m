
% laurene.donati@epfl.ch

% This script runs the ...

clear all; clc; close all;

% Set default paths 
indComp = 4;
switch indComp
    case 1 % Personal lab computer (Mac)
        mainPath = '/Users/laurenedonati/Dropbox/multires-cryo';
    case 2 % Personal laptop (Mac)
        mainPath = '/Users/laurene/Dropbox/multires-cryo';
    case 3 % Masih's laptop (Mac)
        mainPath = '/Users/masih.nilchian/Dropbox/multires-cryo';     
    case 4 % Lab Dell laptop (Linux)
        mainPath = '/home/big/Desktop/multires-cryo';  
end 
addpath(genpath(mainPath)); cd([mainPath '/functions']);


%% INPUTS 

% Paths to data 
volName       = 'ribosome';
versionTest   = 2; 

% Reconstruction parameters
scale         = 1;
alphaRec      = 1e2; %1e-7;
nbItADMM      = 30;
nbItCG        = 10; % 5

% Booleans 
usePositivity = false; 
DISPFIG       = true; 
PRINT         = true; 
SAVEDATA      = false; 

% Initializations 
switch versionTest
    case 1
        pathXmd       = [mainPath '/data/volumes/' volName '/images.xmd'];
        dirStk        = [mainPath '/data/volumes/' volName '/'];
        pathVol       = [mainPath '/data/volumes/' volName '/' volName '.map'];
    case 2
        projectPath   = '/home/big/ScipionUserData/projects/ribosome/';
        pathXmd       = [projectPath 'Runs/000185_XmippProtReconstructHighRes/angles.xmd'];
        dirStk        = projectPath; 
        pathVol       = [mainPath '/data/volumes/' volName '/' volName '.map'];
        symLR         = '/home/big/LR.txt';
end

%% RUN RECONSTRUCTION
[Htb, kernel, rec] = reconstruct_multires_ADMM_xmipp_v2(pathXmd, scale, alphaRec, nbItADMM, nbItCG, usePositivity, symLR); 
%[Htb, kernel, rec] = reconstruct_multires_ADMM_xmipp(pathXmd, dirStk, scale, alphaRec, nbItADMM, nbItCG); 


%% POST-PROCESSING 

% Display results 
if DISPFIG
    vol = xmipp_read(pathVol);
    display_result_multires_rec(vol,rec,size(rec,1));
end 

% Print data to desktop (to visualize in Chimera)
if PRINT
    filePrintRec = '/home/big/Desktop/rec.tif';
    mat2tiffstack(filePrintRec, rec);
    filePrintVol = '/home/big/Desktop/vol.tif';
    mat2tiffstack(filePrintVol, vol);
end


% Save results 
if SAVEDATA
    save([pathXmd 'rec_' volName '_020718_scale' num2str(scale) '.mat']);
end 
 