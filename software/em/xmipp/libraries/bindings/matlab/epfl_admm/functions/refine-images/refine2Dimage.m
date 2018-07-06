
% laurene.donati 
% This function refines 2D image with 1/tau more samples.
% Keeps the middle in the middle 

function refIm = refine2Dimage(im, tau)

% Image sizes 
Nx     = size(im,1); 
Ny     = size(im,2); 

% Original grid 
xVal   = 1:1:Nx; 
yVal   = 1:1:Ny;
[X,Y]  = meshgrid(xVal,yVal);

% Refined grid 
xVal_ref   = 1:tau:Nx; 
yVal_ref   = 1:tau:Ny;
[X_ref,Y_ref]  = meshgrid(xVal_ref,yVal_ref);

% Interpolate image on refined grid
% We can refine here after 
refIm = interp2(X,Y,im,X_ref,Y_ref, 'spline'); 

end 