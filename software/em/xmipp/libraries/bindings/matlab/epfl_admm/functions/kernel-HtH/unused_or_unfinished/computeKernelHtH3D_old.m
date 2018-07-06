%======================================================================== %
% kernel = computeHTH_KB_ConventionalCT_3D_r1r2_clean(param, r1, r2, scale, k1,k2,k3)
%======================================================================== %
% DEFINITION
% This function computes the kernel of HTH using the convolution of the 
% projection of each basis. It performs the computation for Kaiser-Bessel 
% window functions. It computes the convolution using FFT and IFFT, and
% then interpolate it on the grid <k,\theta>.
%======================================================================== %
% INPUTS 
% param.a        support of the KBWF
% param.alpha    window taper parameter of the KBWF
% param.m        smoothness parameter of the KBWF
% r1             def
% r2             def
% scale          scaling parameter of the basis function 
% k1             size of volume along x-dimension
% k2             size of volume along y-dimension
% k3             size of volume along z-dimension
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

function kernel = computeKernelHtH3D_old(param, r1, r2, scale, k1,k2,k3)

% Set kbwf parameters 
a        =   param.a;
alpha    =   param.alpha;
m        =   param.m;

% Compute the 2D projection of one KBWF basis function (put as sub-function)
nberSample = 1001;
x_value    =   linspace(-2*a,2*a,nberSample);
[X,Y]      =   meshgrid(x_value,x_value)    ;
s          =   sqrt(X.^2+Y.^2)              ;
X_KBWF     =   zeros(length(x_value),length(x_value));  
z          =   alpha*sqrt(1-(abs(s)/a).^2);
X_KBWF(abs(s)<a) =  (a/besseli(m,alpha))*sqrt(2*pi/alpha)*...
    ((z(abs(s)<a)/alpha).^(m+1/2)).*...
    besseli(m+1/2,z(abs(s)<a));

% Compute the 2D autocorrelation in Fourier 
yConv    =  ifftshift(ifft2(fft2(X_KBWF).*fft2(X_KBWF)))*16*a^2/((nberSample-1)^2);

% Perform antialiasing
if (param.aliasingFilter)
    yConv = doAntiAliasingFilter2D(yConv, scale*4*a/(nberSample-1));
    
end

% Look-up table with half of 1D line (as auto-correlation is isotropic function)
LT       =   yConv((nberSample+1)/2,(nberSample+1)/2:end);
LTindex  =   x_value((nberSample+1)/2:end);

% Initialize kernel 
kernel           =  zeros(2*k1-1,2*k2-1,2*k3-1);
[K1,K2,K3]       =  meshgrid((-k1+1):(k1-1),(-k2+1):(k2-1),(-k3+1):(k3-1));

% Fill the kernel
for i = 1 : size(r1,2)
    kInerr1          =   K1*r1(1,i)+K2*r1(2,i)+K3*r1(3,i); % inner product <k,r1>
    kInerr2          =   K1*r2(1,i)+K2*r2(2,i)+K3*r2(3,i); % inner product <k,r2>
    KInerAngle       =   sqrt(kInerr1.^2+kInerr2.^2);      % distance from center 
    kernel(KInerAngle<(2*a))  =   kernel(KInerAngle<(2*a))+scale^4*interp1(LTindex,LT,KInerAngle(KInerAngle<(2*a)));
end

end 


