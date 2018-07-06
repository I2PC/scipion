%======================================================================== %
% kernel = computeKernelHtH3D(param, stepSize, r1, r2, scale, sizeRec, ctf)
%======================================================================== %
% DEFINITION
% This function computes the kernel of HTH using the convolution of the
% projection of each basis. It performs the computation for Kaiser-Bessel
% window functions. It computes the convolution using FFT and IFFT, and
% then interpolate it on the grid <k,\theta>.
%======================================================================== %
% INPUTS
% param.a        support of the KBWF
% r1             def
% r2             def
% scale          scaling parameter of the basis function
% k1             size of volume along x-dimension
% k2             size of volume along y-dimension
% k3             size of volume along z-dimension
% varargin       ctf volume (optional)
%======================================================================== %
% OUTPUTS
%  kernel        kernel of HtH
%======================================================================== %
% COPYRIGHT
% Authors:   masih.nilchian@epfl.ch, laurene.donati@epfl.ch
% Date:      June 2018
% Copyright: Biomedical Imaging Group, EPFL
% Citation:  Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.
%======================================================================== %
% WARNING: For the moment, we use the support of the proj of kbwf to
% determine the support of the convolution with the ctf
% WARNING: For the moment, we use the support of the proj of kbwf to
% determine the support of the convolution with the ctf
% Note: Only for isotropic ctf and unique for all projections 
% WARNING : Only work with CTF at scale 1


function kernel = computeKernelHtH3D(param, stepSize, r1, r2, scale, sizeRec, ctf)

disp('-- Computing HtH kernel.')
tic 

% Volume sizes
k1 = sizeRec(1); 
k2 = sizeRec(2); 
k3 = sizeRec(3);

% Support in the spatial domain 
if ~ctf.use 
    suppParam  = 2; 
else 
%     suppParam  = 4; % related to ctf support  
end     

% Compute refined vector of samples  
supp        = suppParam*param.a; 
xVal_acKbwf = -supp:stepSize:supp;
nberSample  = 2*supp/stepSize+1; 

% Compute the refined 2D projection of a KBWF basis function
projKbwfLT2D =  KaiserBesselProjection2D(param, xVal_acKbwf);

if ~ctf.use
    % --- Compute (via Fourier) the 2D autocorrelation of the 2D proj of KBWF
    acKbwfLT2D    =  ifftshift(ifft2(fft2(projKbwfLT2D).*fft2(projKbwfLT2D)))*(stepSize^2);
    % See paper, equ. (16) 
    expVal        = 4; 
elseif ctf.use && scale==1
%     % --- Compute (via Fourier) the 2D autocorrelation of the 2D projection of KBWF and the ctf.   
%     % Compute the autocorrelation (in Fourier) of the proj of KBWF
%     acProjF    = fft2(projKbwfLT2D).*fft2(projKbwfLT2D);
%     % Compute the CTF on a fine grid 
%     fineCTF     = create2Dctf((1/stepSize)+1, nberSample, ctf);
%     % Compute the autocorrelation (in Fourier) of the ctf
%     acCtfF     = fineCTF.*fineCTF;
%     % Perform the point-wise multiplication and inverse FFT
%     acKbwfLT2D      = real(ifftshift(ifft2(acProjF.*acCtfF)))*(stepSize^6);
%     % See paper, equ. (16) 
%     expVal     = 6; 
end

% Look-up table with half of 1D line (as auto-correlation is isotropic function)
centerAC   = (nberSample+1)/2; 
y_LT1D_ac  =  acKbwfLT2D(centerAC,centerAC:end);
x_LT1D_ac  =  xVal_acKbwf(centerAC:end);

% Initialize kernel and meshgrid 
kernel      =  zeros(2*k1-1,2*k2-1,2*k3-1);
[K1,K2,K3]  =  meshgrid((-k1+1):(k1-1),(-k2+1):(k2-1),(-k3+1):(k3-1));

% Fill the kernel inside the support 
for i = 1 : size(r1,2)
    kInner1          =   K1*r1(1,i)+K2*r1(2,i)+K3*r1(3,i); % inner product <K,r1> - relates to y1 
    kInner2          =   K1*r2(1,i)+K2*r2(2,i)+K3*r2(3,i); % inner product <K,r2> - relates to y2 
    KInnerAngle      =   sqrt(kInner1.^2+kInner2.^2);      % distance from center
    kernel(KInnerAngle<(supp))  =  kernel(KInnerAngle<(supp))+scale^expVal*interp1(x_LT1D_ac,y_LT1D_ac,KInnerAngle(KInnerAngle<(supp)));
end

disp(['-- Time to compute kernel: ' num2str(toc) ' sec']);

end


