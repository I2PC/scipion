

% laurene.donati@epfl.ch 
% Test the refinement of a 2D projection image 

clear all
% Set default paths
indComp = 2;
switch indComp
    case 1 % Personal computer (Mac)
        mainPath = '/Users/laurenedonati/Dropbox/multires-cryo';
    case 2 % Personal laptop (Mac)
        mainPath = '/Users/laurene/Dropbox/multires-cryo';
end
addpath(genpath(mainPath));
cd([mainPath '/functions']);

%% Initializations 
volName       = 'betagal';
sizeIndex     = 64;
nberProjs     = 1;
kbwfProj    = struct('a', 2, 'alpha', 10.8, 'm',2);
tauLtRec1D  = 0.01;

coeffPath   = setKBWFCoeffPath(mainPath,volName,sizeIndex);

% Load KBWF coeffs
ckVol    = struct2array(load(coeffPath));

% Set gt size (here only valid for cubic volumes)
sizeVol = [sizeIndex; sizeIndex; sizeIndex];

% Set equi-distributed projection directions
[rot,tilt,psi]  = generateEquidistributedAngles(nberProjs);
[r1,r2]         = setProjectionPlanes(rot,tilt,psi); 

% Compute look-up table with KBWF projection
ltProj1D        = KaiserBesselProjection1D(kbwfProj, tauLtRec1D);

% Compute projection measurements 
disp('-- Compute projections.');
projs    = projection3DMT(ckVol, r1, r2, [1,1,1], [1,1], ltProj1D, kbwfProj.a);



%% Refine projection 
tau    = 0.5;
Nx     = size(projs,1); 
Ny     = size(projs,2); 

xVal   = 1:1:Nx; 
yVal   = 1:1:Ny;
[X,Y]  = meshgrid(xVal,yVal);
%X  = 1:size(projs,1); 
%Y  = 1:size(projs,2); 

xValRef   = 1:tau:Nx; 
yValRef   = 1:tau:Ny;
[Xref,Yref]  = meshgrid(xValRef,yValRef);
%Xq = 1:tau:size(projs,2); 
%Yq = 1:tau:size(projs,2); 
projsRef = interp2(X,Y,projs,Xref,Yref); 
projsRef2 = refine2Dimage(projs,tau); 

% Display
figure, imagesc(projs), colormap gray, title('2D Proj')
figure, imagesc(projsRef), colormap gray, title('Refined 2D Proj')
figure, imagesc(projsRef2), colormap gray, title('Refined 2D Proj -routine ')
