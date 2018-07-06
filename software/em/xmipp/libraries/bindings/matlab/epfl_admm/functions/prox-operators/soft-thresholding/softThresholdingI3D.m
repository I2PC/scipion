% =====================        FUNCTION        ========================== %
% y = softThresholdingI(x,lambda)
% This function performs a soft thresholding on the input:
% y = (norm(x)-lambda)_{+}(x)/norm(x)
% =====================         INPUTS         ========================== %
% x         any real valued input
% lambda    soft thresholding parameter 
% =====================         OUTPUTS        ========================== %
% y         soft-thresholded output (same size as x) 
%======================================================================== % 

%======================================================================== %
% Authors:  masih.nilchian@epfl.ch, laurene.donati@epfl.ch
% Date: June 2018
% Copyright: Biomedical Imaging Group, EPFL
% Citation: Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.
%======================================================================== %

function y = softThresholdingI3D(x,lambda)
x1 = x;
if ndims(x)==2
    x=reshape(x,[round(size(x,1)^(1/3)),round(size(x,1)^(1/3)),round(size(x,1)^(1/3)),3]);
end
    
xNorm          = sqrt(sum(x.^2,4))          ;
xNorm          = repmat(xNorm,[1,1,1,3]);
y              = xNorm                      ;
y(xNorm~=0)    = y(xNorm~=0)-lambda         ;
y(y<0)         = 0                          ;
y(xNorm~=0)    = y(xNorm~=0).*x(xNorm~=0)./xNorm((xNorm~=0)) ;
y              = reshape(y,size(x1))        ;