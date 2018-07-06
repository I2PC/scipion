% laurene.donati@epfl.ch

% Format projections and directions from Xmipp output s 

function [projs, rot, tilt, psi, shiftX, shiftY, flipList] = formatProjsInfoXmipp(pathXmd)


% Read xmd file and go to directory 
%s = xmipp_read_metadata([pathXmd 'images.xmd']); 
s = xmipp_read_metadata(pathXmd); 
%cd(dirStk);

selected = find(s.enabled==1); 
% % TEST FROM HERE FOR FLIPPING
% flipList = s.flip(selected);
% selected = selected(flipList==1); 
% flipList = flipList(flipList==1);
% % TO HERE 

% Get sizes 
nberProjs = length(selected);

%nberProjs   = 15; 
projSize  = size(xmipp_read(s.image{1}));

% Initializations 
projs      = zeros([projSize, nberProjs]);

% Format projection measurements 
ind = 1; 
for p = selected
    projs(:,:,ind) = xmipp_read(s.image{p});
    ind = ind+1; 
end 

% Return projection angles in radians 
rot  = (pi/180)*s.angleRot(selected); 
tilt = (pi/180)*s.angleTilt(selected);
psi  = (pi/180)*s.anglePsi(selected);

% Return shifts information
shiftX = s.shiftX(selected);
shiftY = s.shiftY(selected);

% Return mirror information
flipList = s.flip(selected);

end 