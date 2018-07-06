
% laurene.donati@epfl.ch 
% This function applies multiresolution reconstruction to xmipp projections

function [Htb, kernel, rec] = reconstruct_multires_ADMM_xmipp_v2(pathXmd, scale, alphaRec, nbItADMM, nbItCG, usePositivity, symLR)

% KBWF parameters (for reconstruction)
kbwfRec = struct('a', 4, 'alpha', 19, 'm', 2);

% CTF information
ctf.use = false;

% Process projection measurements and get Htb, r1, r2 
[Htb, r1, r2, sizeVol, sizeRec] = processProjsAndComputeHtbXmipp(pathXmd, kbwfRec, scale, symLR); 

% Compute kernel and define HtH (with kernel)
tauKernel = 0.25; % (To optimize)
kernel = computeKernelHtH3D(kbwfRec,tauKernel, r1, r2, scale, sizeRec, ctf);
HTH    = @(X) applyKernel3D(kernel,X);

% Set regularization (prox operator)
prox = @(X,alphaRec) softThresholdingI3D(X,alphaRec);

% Set projection constraints
if usePositivity 
    projConstr = @(X) positiveProj(X);
else
   projConstr  = @(X) X;
end 

tic 
% Run ADMM reconstruction 
[ckRec,~] = applyADMM(Htb,'ATA',HTH,'verbose', true,'drawflag',false, ...
        'alpha',alphaRec,'lambda',0,'tol',1e-20,'maxOutIter', nbItADMM, ...
    'maxInIter',nbItCG, 'mu',1e1*alphaRec,'doProjection',projConstr,'prox',prox);

% Re-expand reconstructed ck in fine grid
rec = expandScaledKBWFcoeffstoFineGrid_3D(ckRec, scale, kbwfRec, sizeVol);

% Correct for misflipped reconstruction 
%rec = permute(rec, [2 1 3]);

% Display reconstruction time 
disp(['Time to run ADMM and re-expand: ' num2str(toc) ' sec'])

end 