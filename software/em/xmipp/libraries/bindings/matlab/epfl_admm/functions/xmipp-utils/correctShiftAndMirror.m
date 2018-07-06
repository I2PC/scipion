
% laurene.donati@epfl.ch
% This function corrects the shift of the projections through spline interpolation
% The shift information are output from Scipion.
% NOTE: Currently works only for square projection images.

function projs = correctShiftAndMirror(projs,shiftX,shiftY,flipList)


% Padd projections to handle shift at border
maxShift    = ceil(max(abs([shiftX(:); shiftY(:)])));
paddedProjs = padarray(projs,[maxShift,maxShift],0,'both');

% Get size info
Nx  = size(projs,1);
Nxp = size(paddedProjs,1);
P   = size(projs,3);

% Get centered projection meshgrid
cP        = (Nxp-1)/2;
xValP     = -cP:cP;
yValP     = -cP:cP;
[pGx,pGy] = meshgrid(xValP,yValP);


%% INTERPOLATION (SHOULD BE PARALLELISED)
for p = 1:P
    
    % Get interpolation center
    cI = (Nx-1)/2;
    
    % Handle shifts and projection mirrors (TO TEST!!!)
    if ~flipList(p)
        xValI     = (-cI:cI) - shiftX(p);
        yValI     = (-cI:cI) - shiftY(p);
    else
        xValI     = -(-cI:cI) + shiftX(p);
        yValI     = (-cI:cI) - shiftY(p);
    end
    
    % Create interpolation meshgrid 
    [iGx,iGy] = meshgrid(xValI,yValI);
    
    % Interpolate with sinc
    projs(:,:,p) = interp2(pGx,pGy,paddedProjs(:,:,p),iGx,iGy, 'sinc');
    
end

end

% TESTS     
% TEST = true;
%     projsPost(:,:,p) = interp2(pGx,pGy,paddedProjs(:,:,p),iGx,iGy);
%     if TEST && p<11
%         figure, subplot(2,1,1), imagesc(projs(:,:,p)), grid on;  title('Before shift correction');
%         subplot(2,1,2), imagesc(projsPost(:,:,p)), grid on; ...
%             title(['After (shiftX=' num2str( shiftX(p)) ', shiftY=' num2str(shiftY(p)) ')']); 
%     end