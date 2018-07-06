% =====================        FUNCTION        ========================== %
% [r1, r2] = setAngularGeometry(rot,tilt,psi)
%======================================================================== %
% DEFINITION
% This function gives the two first lines of the 3x3 rotation matrix that
% converts the object coordinate system to the projection coordinate system. 
% Note: These first two lines permit to obtain the coordinates of the 
% projection of a 3D object on the projection planes.
% Note: The third line r3 (projection direction) is not needed as the value 
% of a projection along this axis is always zero by definition. 
%======================================================================== %
% INPUTS 
% rot     Angle between the x axis and the N axis [in rad].
% tilt    Angle between the z axis and the Z axis [in rad].
% psi     Angle between the N axis and the X axis. Zero by default. 
%======================================================================== % 
% OUTPUTS 
% r1      Vector that returns the value of a 3D point of the object in the 
%         first projection coordinate axis y1 defined by the Euler angles. 
% r2      Vector that returns the value of a 3D point of the object in the 
%         second projection coordinate axis y2 defined by the Euler angles.  
%======================================================================== %

%=========================================================================
% COPYRIGHT 
% Authors: masih.nilchian@epfl.ch, laurene.donati@epfl.ch
% Date: June 2018
% Copyright: Biomedical Imaging Group, EPFL
% Citation: Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.
%=========================================================================

function [r1,r2] = setProjectionPlanes(rot,tilt,psi,useEulerXmipp)

% Initializations
nberProj  = length(rot);
r1        = zeros(3,nberProj);
r2        = zeros(3,nberProj);

% Computing rotation matrix
for i = 1:nberProj
    if useEulerXmipp
        R = Euler(rot(i), tilt(i), psi(i));
    else
        R = Euler_xmipp(rot(i), tilt(i), psi(i));
    end 
    r1(:,i) = [R(1,1),R(1,2),R(1,3)]';
    r2(:,i) = [R(2,1),R(2,2),R(2,3)]';
end

end