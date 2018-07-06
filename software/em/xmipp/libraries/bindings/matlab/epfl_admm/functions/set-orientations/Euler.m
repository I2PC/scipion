% =====================        FUNCTION        ========================== %
% R = Euler(phi, theta, psi)
% This function rotates the volume axis to match the projection axis.  
% The first projection plane vector (y1) is given along the x direction 
% The second projection plane vector (y2) is given along the y direction 
% The projection direction vector (y3) is given along the z direction
% =====================         INPUTS         ========================== %
% rot     def
% tilt    def
% psi     def
% =====================         OUTPUTS        ========================== %
% R      def
%======================================================================== %

%=========================================================================
% Authors: masih.nilchian@epfl.ch, laurene.donati@epfl.ch
% Date: June 2018
% Copyright: Biomedical Imaging Group, EPFL
% Citation: Fast Multiresolution Reconstruction for Cryo-EM, Donati et al.
%=========================================================================

function R = Euler(phi, theta, psi)

	R = [cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi), ...  % R(1,1)
         cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi), ...  % R(1,2)
         -cos(psi)*sin(theta); ...                            % R(1,3)
		 -sin(psi)*cos(theta)*cos(phi)-cos(psi)*sin(phi), ... % R(2,1)
         -sin(psi)*cos(theta)*sin(phi)+cos(psi)*cos(phi), ... % R(2,2)
         sin(psi)*sin(theta); ...                             % R(2,3)
		 sin(theta)*cos(phi), ...                             % R(3,1)
         sin(theta)*sin(phi), ...                             % R(3,2)
         cos(theta)];                                         % R(3,3)

end