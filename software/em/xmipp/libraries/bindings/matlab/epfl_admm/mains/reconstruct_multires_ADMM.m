
function rec = reconstruct_multires_ADMM(projs, r1, r2, scale, alphaRec, nbItADMM, nbItCG, sizeIndex)

% QU1: Add here conversion from (rot, tilt, psi) to (r1, r2) ? 
% [r1,r2] = setProjectionPlanes(rot,tilt,psi);

% KBWF parameters (for reconstruction)
kbwfRec     = struct('a', 4, 'alpha', 19, 'm', 2);

% CTF information
ctf.use     = false;

% Refinement paramaters for look-up tables
tauHtb      = scale*0.25; 
tauKernel   = 0.01;

% Set reconstruction size (based on scale)
sizeVol    = [sizeIndex; sizeIndex; sizeIndex];
sizeRec    = floor((sizeVol-1)/scale+1);

% Combut Htb (fast version)
Htb        = computeHtbfast(projs, kbwfRec, tauHtb, r1, r2, scale, sizeRec, 0);

% Compute kernel and define HtH (with kernel)
kernel     = computeKernelHtH3D(kbwfRec,tauKernel, r1, r2, scale, sizeRec, ctf);
HTH        = @(X) applyKernel3D(kernel,X);

% Set regularization (prox operator)
prox       = @(X,alphaRec) softThresholdingI3D(X,alphaRec);

% Set projection constraints
projConstr = @(X) positiveProj(X);

% Run ADMM reconstruction 
tic
[ckRec,~]  = applyADMM(Htb,'ATA',HTH,'verbose', true,'drawflag',false, ...
        'alpha',alphaRec,'lambda',0,'tol',1e-20,'maxOutIter', nbItADMM, ...
    'maxInIter',nbItCG, 'mu',1e1*alphaRec,'doProjection',projConstr,'prox',prox);

% Re-expand reconstructed ck in fine grid
rec        = expandScaledKBWFcoeffstoFineGrid_3D(ckRec, scale, kbwfRec, sizeVol);

% Display reconstruction time 
disp(['Time to run ADMM and re-expand: ' num2str(toc) ' sec'])

end 