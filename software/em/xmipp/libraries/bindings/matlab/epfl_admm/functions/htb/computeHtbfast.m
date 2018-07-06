%======================================================================== %
% Htb = computeHtbfast(projs, kbwf_p, tau, r1, r2, scale, sizeRec, padding)
%======================================================================== %
% DEFINITION
% This function computes the product HTb based on equation (17).
% It performs the computation for scaled Kaiser-Bessel window functions.
% It computes the convolutions using FFT and IFFT, and then interpolates
% it on the grid <k,\theta>.
%======================================================================== %
% INPUTS
% projs          Set of 2D projection measurements
% kbwf_p        Structure containing the KBWF parameters
% tau 
% r1             def
% r2             def
% scale          Scaling parameter of the basis function
% sizeRec        Size of volume along [x, y, z]
% padding        
%======================================================================== %
% OUTPUTS
% Htb            Adjoint of H applied to the projections
%======================================================================== %

%======================================================================== %
% COPYRIGHT
% Authors:   laurene.donati@epfl.ch, masih.nilchian@epfl.ch
% Date:      June 2018
% Copyright: Biomedical Imaging Group, EPFL
% Citation:  Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.
%======================================================================== %

% Note, we do not consider ctf yet. For now, only works for volumes of 
% even sizes AND cubic.

function Htb = computeHtbfast(projs, kbwf_p, tau, r1, r2, scale, sizeRec, padding)

disp('-- Computing Htb.');
tic 

% Get size information
k1 = sizeRec(1); 
k2 = sizeRec(2); 
k3 = sizeRec(3);
nberProjs = size(projs,3);

% Compute FT of refined 2D KBWF projection 
kbwf   = KaiserBesselProjection2D_RefinedScaled(kbwf_p,tau,scale,size(projs));
kbwf   = padarray(kbwf,round(size(kbwf)*padding),0,'both'); 
fkbwf  = fft2(kbwf);

% Initialize Htb and its meshgrid
Htb         = zeros(k1,k2,k3);
centerHtb   = (k1-1)/2;
xValHtb     = scale*(-centerHtb:centerHtb); 
[K1,K2,K3]  = meshgrid(xValHtb, xValHtb, xValHtb);


%% For each particle p (TO BE PARALLELISED) :

for p=1:nberProjs
    
    % Compute measurement on refined grid
    proj_p      = refine2Dimage(projs(:,:,p), tau);
    proj_p      = padarray(proj_p,round(size(proj_p)*padding),0,'both');
    
    % Compute the convolution as point-wise multiplication in Fourier
    conv_p      = ifftshift(ifft2(fkbwf.*fft2(proj_p)))*(tau^2);
    
    % Set meshgrid in projection domain 
    centerConv  = (size(conv_p,1)-1)/2;
    xValConv    = tau*(-centerConv:centerConv);
    [Y1,Y2]     = meshgrid(xValConv, xValConv);
    
    % Fill Htb
    kInner1     = K1*r1(1,p)+K2*r1(2,p)+K3*r1(3,p); % inner product <K,r1> - relates to y1
    kInner2     = K1*r2(1,p)+K2*r2(2,p)+K3*r2(3,p); % inner product <K,r2> - relates to y2
    contrib_p   = interp2(Y1, Y2, conv_p, kInner1(:), kInner2(:), 'spline');
    Htb         = Htb + reshape(contrib_p,size(Htb));
end

disp(['-- Time to compute Htb: ' num2str(toc) ' sec']);

end
