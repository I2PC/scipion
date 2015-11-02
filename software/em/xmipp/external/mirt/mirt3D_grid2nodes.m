% MIRT3D_GRID2NODES  transforms the dense gradient to the
% gradient at each node (B-spline control point).
%
% Input 
% ddx, ddy, ddz - the gradient of similarity measure at each image voxel
% in x, y and z direction respectively
%
% F - is a precomputed matrix of B-spline basis function coefficients, see.
% mirt3D_F.m file
%
% okno - is a spacing width between the B-spline control points
% mg, ng, kg - size of B-spline control points in x,y and z directions
%
% Output
% Gr - a 4D array of the gradient of the similarity measure at B-spline
% control point positions. The first 3 dimenstion (size = mg x ng x kg) is the organized control
% points. The 4th dimesion (size 3) is the index of x,y or z component.

% Copyright (C) 2007-2010 Andriy Myronenko (myron@csee.ogi.edu)
% also see http://www.bme.ogi.edu/~myron/matlab/MIRT/
% This file is a part of Medical Image Registration Toolbox (MIRT)


function Gr=mirt3D_grid2nodes(ddx, ddy, ddz, F, okno, mg, ng, kg)

% greate the 3D matrices wich includes
% the gradient of the similarity measure at
% each of the B-spline control points
grx=zeros(mg,ng,kg);
gry=zeros(mg,ng,kg);
grz=zeros(mg,ng,kg);

% cycle over all 4x4x4 patches of B-spline control points
for i=1:mg-3,
    for j=1:ng-3,
        for k=1:kg-3,
            
        % define the indices of the voxels corresponding
        % to the given 4x4x4 patch    
        in1=(i-1)*okno+1:i*okno;
        in2=(j-1)*okno+1:j*okno;
        in3=(k-1)*okno+1:k*okno;
        
        % extract the voxel-wise gradient (in x direction) of the similarity measure
        % that correspond to the given 4x4x4 patch of B-spline cotrol
        % points. 
        tmp=ddx(in1,in2,in3);
        
        % multiply the voxel-wise gradient by the transpose matrix of
        % precomputed B-spline basis functions and accumulate into node-wise (B-spline)
        % gradient. Accumulation is because some B-spline control
        % points are shared between the patches
        grx(i:i+3,j:j+3,k:k+3,1)=grx(i:i+3,j:j+3,k:k+3,1)+reshape(F'*tmp(:),[4 4 4]);
      
        % do the same thing for y and z coordinates
        tmp=ddy(in1,in2,in3);
        gry(i:i+3,j:j+3,k:k+3,1)=gry(i:i+3,j:j+3,k:k+3,1)+reshape(F'*tmp(:),[4 4 4]);
      
        tmp=ddz(in1,in2,in3);
        grz(i:i+3,j:j+3,k:k+3,1)=grz(i:i+3,j:j+3,k:k+3,1)+reshape(F'*tmp(:),[4 4 4]);
         
       end 
    end
end

% concatinate into a single 4D array
Gr=cat(4, grx, gry,grz);