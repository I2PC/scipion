% ========== HERE TEST FOR THE KERNEL WITH CTF ==========
TEST_kernel = true;
if TEST_kernel
    testIm   = zeros(sizeIndex,sizeIndex,sizeIndex);
    center = floor(sizeIndex/2);
    testIm(center, center, center) =1;
    projs    = projection3DMT(testIm, r1, r2, [1,1,1], [1,1], ltProj, kbwfProj.a);
    if ctf.use == true
        ctf.image     = create2Dctf(size(projs), size(projs), ctf);
    end
    % === TEST H^T(ctf(ctf(H(volCk)))) ===
    % Up to here, projs is H(volCk)
    projs     = applyCTFAdjoint(projs,ctf);
    % Up to here, projs is ctf(H(volCk))
    projs     = convolveProjsWithCTF(projs,ctf);
    % Up to here, projs is ctf(ctf(H(volCk)))
    sizeRec   = floor((sizeVol-1)/scale+1);
    ltRec1D   = KaiserBesselProjection1D(kbwfRec, tauLtRec1D);
    resHtH   = projectionAdjoint3DMT(projs, sizeRec, r1, r2, ...
        scale*[1,1,1], [1,1], ltRec1D, kbwfRec.a);
    % Here, resHtH is H^T(ctf(ctf(H(volCk))))
    if ctf.use
        titleHtH = 'HtHc (with CTF)';
    else
        titleHtH = 'HtHc (no CTF)';
    end
    figure, imagesc(resHtH(:,:,sizeIndex/2)), colormap gray, title(titleHtH);
    
    
    % === TEST kernel*ckVol ===
    kernel    = computeKernelHtH3D_v2(kbwfRec,tauLtRec2D, r1, r2, scale, sizeRec, ctf);
    resKern = applyKernel3D(kernel,testIm);
    if ctf.use
        titleKernel = 'r*c (with CTF)';
    else
        titleKernel = 'r*c (no CTF)';
    end
    figure, imagesc(resKern(:,:,sizeIndex/2)), colormap gray, title(titleKernel);
    
end

% ================================================================================

