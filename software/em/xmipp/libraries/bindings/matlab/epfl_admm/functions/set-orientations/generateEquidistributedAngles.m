% =====================        FUNCTION        ========================== %
% [rot, tilt, psi] = generateEquidistribRandomProjAngles_clean(nberProj)
%======================================================================== %
% DEFINITION
% This function randomly generates a number nberProj of  equidistributed 
% 3D projection (Euler) angles for simulation purposes. 
% Based on [Deserno, Markus. "How to generate equidistributed ...", 2004]
% Self-note: In our notations, his phi is our rot, his theta is our tilt. 
%======================================================================== %
% INPUTS
% nberProj    Number of total projection directions.         
%======================================================================== %
% OUTPUTS
% rot         Angle between the x axis and the N axis [rad].
% tilt        Angle between the z axis and the Z axis [rad].
% psi         Angle between the N axis and the X axis. Zero by default. 
%======================================================================== % 

%======================================================================== %
% COPYRIGHT
% Authors: laurene.donati@epfl.ch, masih.nilchian@epfl.ch
% Date: June 2018
% Copyright: Biomedical Imaging Group, EPFL
% Citation: Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.
%=========================================================================

function [rot, tilt, psi] = generateEquidistributedAngles(nberProj)

 % Set rot angles 
 z    = (rand(1,nberProj)*2)-1;
 rot  = rand(1,nberProj)*2*pi;
 
 % Set tilt angles
 x    = sqrt(1-(z.^2)).*cos(rot);
 y    = sqrt(1-(z.^2)).*sin(rot);
 tilt = atan(sqrt(y.^2+x.^2)./z);
 tilt(tilt<0) = tilt(tilt<0)+pi; 
 
 % Set psi angles to zero 
 psi = zeros(size(rot)); 

end 