%======================================================================== %
% Htb = computeHtbfast(projs, kbwfRec, ltRec1D, r1, r2, scale, sizeRec)
%======================================================================== %
% DEFINITION
% This function computes the product HTb based on equation (17).
% It performs the computation for scaled Kaiser-Bessel window functions.
% It computes the convolutions using FFT and IFFT, and then interpolates
% it on the grid <k,\theta>.
%======================================================================== %
% INPUTS
% projs          Set of 2D projection measurements
% kbwfRec        Structure containing the KBWF parameters
% ltRec1D        1D projection of the KBWF  ???????
% r1             def
% r2             def
% scale          Scaling parameter of the basis function
% sizeRec        Size of volume along [x, y, z]
% ctf            (to be added)
%======================================================================== %
% OUTPUTS
% Htb            Adjoint of H applied to the projections.
%======================================================================== %
% COPYRIGHT
% Authors:   laurene.donati@epfl.ch, masih.nilchian@epfl.ch
% Date:      June 2018
% Copyright: Biomedical Imaging Group, EPFL
% Citation:  Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.
%======================================================================== %

% Warning, we do not consider scaling and ctf yet
% Warning: For now, only works for volumes of even sizes AND cubic
% Question: Keep padding before FFT ? And cropping after?

function Htb = computeHtbfast_v2(projs, kbwf_p, tau, r1, r2, scale, sizeRec, paddRatio)

% Get size information
k1 = sizeRec(1); 
k2 = sizeRec(2); 
k3 = sizeRec(3);
px     = size(projs,1); 
py     = size(projs,2); 
nProjs = size(projs,3);

% Compute centered grid (for measurements) 
[~, ~, xG, yG]                    = set2Dgrid([px; py], 1, 1); 

% Compute~ centered grid after scaling and refinement (for interp + kwbf)
[~, ~, xGRefScaled, yGRefScaled]  = set2Dgrid([px; py], tau, scale); 

% Compute~ centered grid after scaling and refinement (for interp + kwbf)
[xValRef, ~, Y1, Y2]  = set2Dgrid([px; py], tau, 1); 

% Compute FT of padded 2D proj of refined KBWF
kbwf   = KaiserBesselProjection2D(kbwf_p, xValRef);
kbwf   = padarray(kbwf,round(size(kbwf)*paddRatio),0,'both'); %reduce padding !!
fkbwf  = fft2(kbwf);

% Initialize Htb and its meshgrid
Htb         = zeros(k1,k2,k3);
centerHtb   = (k1-1)/2;
[K1,K2,K3]  = meshgrid(-centerHtb:centerHtb,-centerHtb:centerHtb,-centerHtb:centerHtb);


%% For each particle p (note: could be parallelized) :

for p=1:nProjs
    % Interpolate measurement on refined/scaled grid and padd
    proj_p   = interp2(xG,yG,projs(:,:,p),xGRefScaled,yGRefScaled, 'spline'); 
    proj_p   = padarray(proj_p,round(size(proj_p)*paddRatio),0,'both'); %reduce padding
    
    % Compute the convolution as point-wise multiplication in Fourier
    conv_p      = ifftshift(ifft2(fkbwf.*fft2(proj_p)))*((tau*scale)^2);
    %centerConv  = (size(conv_p,1)-1)/2;
    %[Y1,Y2]     = meshgrid(tau*(-centerConv:centerConv),tau*(-centerConv:centerConv));
    % Then crop FFT ?
    
    % Fill Htb
    kInner1   = K1*r1(1,p)+K2*r1(2,p)+K3*r1(3,p); % inner product <K,r1> - relates to y1
    kInner2   = K1*r2(1,p)+K2*r2(2,p)+K3*r2(3,p); % inner product <K,r2> - relates to y2
    contrib_p = interp2(Y1, Y2, conv_p, kInner1(:), kInner2(:), 'spline');
    Htb       = Htb + reshape(contrib_p,size(Htb));
end
% Take into consideration scale effect
Htb = scale* Htb; 
end

% verbose = false; 
% if p==1 && verbose; tic; end
%     if p==1 && verbose; disp(['--- for a particle, time for refining the 2D image: ' num2str(toc) ' sec']); end
%     if p==1 && verbose; disp(['--- for a particle, time for computing the convolution in Fourier: ' num2str(toc) ' sec']); end
% 
%     if p==1 && verbose; disp(['--- for a particle, time for picking and filling the kernel: ' num2str(toc) ' sec']); end

% sizeProj     = ; 
% sizeDSProj   = floor((sizeProj-1)/scale+1);
% centerDSProj = (sizeDSProj(1)+1)/2;
% xValRef      = -centerDSProj:tau:centerDSProj;
% centerProj = (sizeProj(1)-1)/2;
% xValRef      = -centerProj:1:centerProj;

