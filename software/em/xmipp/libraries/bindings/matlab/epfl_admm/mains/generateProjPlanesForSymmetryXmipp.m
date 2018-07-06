
% laurene.donati@epfl.ch

function [r1, r2] = generateProjPlanesForSymmetryXmipp(s,nberProj,nSymm,Rsymm) 

% Initializes r1 and r2
r1  = zeros(3,nberProj*nSymm);
r2  = zeros(3,nberProj*nSymm);

for p = 1:nberProj
    % Read the directions of p-th particle
    rot    = (pi/180)*s.angleRot(p);
    tilt   = (pi/180)*s.angleTilt(p);
    psi    = (pi/180)*s.anglePsi(p);
    % Generate default projection plane
    R0     = Euler(rot, tilt, psi);
    % Set starter index and counter 
    numP    = ((p-1)*nSymm)+1;
    counter = 0;
    for i = 1:nSymm
        % Compute symmetrized matrix
        R = R0*Rsymm(:,:,i)';
        % Set projection planes
        r1(:, numP+counter) = [R(1,1),R(1,2),R(1,3)]';
        r2(:, numP+counter) = [R(2,1),R(2,2),R(2,3)]';
        % Update counter
        counter = counter+1;
    end
end

end 