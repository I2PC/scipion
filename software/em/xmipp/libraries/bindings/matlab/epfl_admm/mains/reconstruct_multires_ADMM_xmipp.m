
% laurene.donati@epfl.ch 
% This function applies multiresolution reconstruction to xmipp projections

function [Htb, kernel, rec] = reconstruct_multires_ADMM_xmipp(pathXmd, scale, alphaRec, nbItADMM, nbItCG)
%function [Htb, kernel, rec] = reconstruct_multires_ADMM_xmipp(pathXmd, dirStk, scale, alphaRec, nbItADMM, nbItCG)

DISPFIG = false; 

% Format Xmipp projection and directions from xmd file 
% ASPECTS TO PARALLELISE 
[projs, rot, tilt, psi, shiftX, shiftY, flipList] = formatProjsInfoXmipp(pathXmd);

% Correct shifts 
% ASPECTS TO PARALLELISE 
projs = correctShiftAndMirror(projs,shiftX,shiftY,flipList);

% Generate projection planes 
useEulerXmipp = false; 
[r1,r2] = setProjectionPlanes(rot,tilt,psi,useEulerXmipp);

% KBWF parameters (for reconstruction)
kbwfRec = struct('a', 4, 'alpha', 19, 'm', 2);

% CTF information
ctf.use = false;

% Refinement paramaters for look-up tables (to optimize)
tauHtb    = scale*0.25; 
tauKernel = 0.25;

% Set reconstruction size (based on scale)
sizeIndex = size(projs,1);
sizeVol   = [sizeIndex; sizeIndex; sizeIndex];
sizeRec   = floor((sizeVol-1)/scale+1);

% Combut Htb (fast version)
% ASPECTS TO PARALLELISE 
Htb = computeHtbfast(projs, kbwfRec, tauHtb, r1, r2, scale, sizeRec, 0.05); % padding ratio for proteasome: 0.05
if DISPFIG; figure, imagesc(Htb(:,:,floor(sizeRec(1)/2))), title('Htb Fast'); end 

% Compute kernel and define HtH (with kernel)
kernel = computeKernelHtH3D(kbwfRec,tauKernel, r1, r2, scale, sizeRec, ctf);
HTH    = @(X) applyKernel3D(kernel,X);
if DISPFIG; figure, imagesc(kernel(:,:,sizeRec(1))), title('kernel'); end 

% Set regularization (prox operator)
prox = @(X,alphaRec) softThresholdingI3D(X,alphaRec);

% Set projection constraints
projConstr = @(X) positiveProj(X);

% Run ADMM reconstruction 
tic
[ckRec,~] = applyADMM(Htb,'ATA',HTH,'verbose', true,'drawflag',false, ...
        'alpha',alphaRec,'lambda',0,'tol',1e-20,'maxOutIter', nbItADMM, ...
    'maxInIter',nbItCG, 'mu',1e1*alphaRec,'doProjection',projConstr,'prox',prox);

% Re-expand reconstructed ck in fine grid
rec = expandScaledKBWFcoeffstoFineGrid_3D(ckRec, scale, kbwfRec, sizeVol);

% Correct for misflipped reconstruction 
rec = permute(rec, [2 1 3]);

% Display reconstruction time 
disp(['Time to run ADMM and re-expand: ' num2str(toc) ' sec'])

end 