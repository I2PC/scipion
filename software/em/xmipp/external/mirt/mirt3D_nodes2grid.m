% MIRT3D_NODES2GRID  computes the dense transformation (positions of all image voxels)
% from the current positions of B-spline control points
%
% Input
% X - 5D array of B-spline control point positions. The first 3 dimensions
% include the coordinates of B-spline control points in a particular
% volume. The 4th dimension is of size 3, it indicates wheather it is the
% X, Y or Z coordinate. The 5th dimension is time (volume number).
%
% F - is a precomputed matrix of B-spline basis function coefficients, see.
% mirt3D_F.m file
%
% okno - is a spacing width between the B-spline control points
%
% Output
% Xx,Xy,Xz - 3D matrices of positions of all image voxels computed from the
% corresponding positions of B-spline control points

% Copyright (C) 2007-2010 Andriy Myronenko (myron@csee.ogi.edu)
% also see http://www.bme.ogi.edu/~myron/matlab/MIRT/

function [Xx,Xy,Xz]=mirt3D_nodes2grid(X, F, okno)

[mg,ng,kg, tmp]=size(X);

% Compute the dense transformation in Xx,Xy,Xz
% Go over all 4x4x4 patches of B-spline control points
for i=1:mg-3,
    for j=1:ng-3,
      for k=1:kg-3,  
        
        % define the indices of the voxels corresponding
        % to the given 4x4x4 patch
        in1=(i-1)*okno+1:i*okno;
        in2=(j-1)*okno+1:j*okno;
        in3=(k-1)*okno+1:k*okno;  
        
        % take the X coordinates of the current 4x4x4 patch of B-spline
        % control points, rearrange in vector and multiply by the matrix
        % of B-spline basis functions (F) to get the dense coordinates of 
        % the voxels within the given patch
        tmp=X(i:i+3,j:j+3,k:k+3,1);
        Xx(in1,in2,in3)=reshape(F*tmp(:),[okno okno okno]);

        % repeat for Y coordinates of the patch
        tmp=X(i:i+3,j:j+3,k:k+3,2);
        Xy(in1,in2,in3)=reshape(F*tmp(:),[okno okno okno]);

        % repeat for Z coordinates of the patch
        tmp=X(i:i+3,j:j+3,k:k+3,3);
        Xz(in1,in2,in3)=reshape(F*tmp(:),[okno okno okno]);
      end
    end
end
