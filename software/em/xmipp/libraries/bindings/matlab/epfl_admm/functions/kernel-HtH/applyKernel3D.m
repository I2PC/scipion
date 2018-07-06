% =====================        FUNCTION        ========================== %
% y = applyKernel3D(kernel,x)
% This function returns the convolution of the input with the kernel
% =====================         INPUTS         ========================== %
% x         3D real valued input
% kernel    HtH kernel  
% =====================         OUTPUTS        ========================== %
% y         convolution of the input with the kernel
%======================================================================== % 

%=========================================================================
% Authors:  masih.nilchian@epfl.ch, laurene.donati@epfl.ch
% Date: June 2018
% Copyright: Biomedical Imaging Group, EPFL
% Citation: Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.
%=========================================================================

function y = applyKernel3D(kernel,x)

% get image size
k1  =  size(x,1);
k2  =  size(x,2);
k3  =  size(x,3);

% apply filter using fft
y   = real(ifftn(fftn(x,[2*k1-1,2*k2-1,2*k3-1]).*fftn(ifftshift(kernel))));

% crop the region of interest
y   = y(1:k1,1:k2,1:k3);

end 
