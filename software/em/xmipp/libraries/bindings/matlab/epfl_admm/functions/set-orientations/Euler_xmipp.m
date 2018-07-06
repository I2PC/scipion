% =====================        FUNCTION        ========================== %
% R = Euler_xmipp(phi, theta, psi)
% This function is an adaptation of Euler.m for handling xmipp data
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

function R = Euler_xmipp(phi, theta, psi)

	R = [cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi), ...
         cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi), ...
         -cos(psi)*sin(theta); ...
		 -sin(psi)*cos(theta)*sin(phi)+cos(psi)*cos(phi), ...
         -sin(psi)*cos(theta)*cos(phi)-cos(psi)*sin(phi), ...
         sin(psi)*sin(theta); ...
		 sin(theta)*cos(phi), ...
         sin(theta)*sin(phi), ...
         cos(theta)];
end