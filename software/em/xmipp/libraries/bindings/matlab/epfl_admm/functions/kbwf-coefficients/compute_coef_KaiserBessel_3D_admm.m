%%   function  Ck = compute_coef_KaiserBessel_3D_admm(volume,alpha,a,m, nberItCoeff)
%   computes 3d coef of the input volume
% SELF NOTE: I HAVE ADDED THE SCALE AS AN OPTION !!!!

%%
function  ck = compute_coef_KaiserBessel_3D_admm(vol,scale, alpha,a,m, nberItCoeff, savePath)

% The way that we fill the coarse grid to fine grid is as follows:
%    fix the boundaries, and refine the gaps in between
%
 H = @(x) Reconstruct_KaiserBessel_3D(x,scale,alpha,a,m,size(x));
    HT= @(x) Reconstruct_KaiserBessel_3D(x,scale,alpha,a,m,size(x));
    HTH= @(x) HT(H(x))+1e-15*x;

ck = applyADMM(HT(vol),'ATA',HTH,'verbose',true,'maxOutIter',nberItCoeff);

save(savePath, 'ck');

% Test the nIter for KBWF coefficient computation
disp('Maybe you should test the coefficients later.')
KBWFcoeffTest = false;
if KBWFcoeffTest
    kbwfTest = Reconstruct_KaiserBessel_3D_center(volCk, scale, ...
        paramProj.alpha,paramProj.a,paramProj.m,[imSize, imSize, imSize]);
end

end
