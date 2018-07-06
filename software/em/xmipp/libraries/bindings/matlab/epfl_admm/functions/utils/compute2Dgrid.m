
% laurene.donati@epfl.ch 
% refine a 2D grid based on refinement parameter tau

function [Gx, Gy, x_val, y_val] = compute2Dgrid(Nx,Ny,tau)

% Compute grid
x_val = -Nx/(2*tau):1:Nx/(2*tau);
y_val = -Nx/(2*tau):1:Ny/(2*tau);
[Gx, Gy]  = meshgrid(x_val,y_val);

end
