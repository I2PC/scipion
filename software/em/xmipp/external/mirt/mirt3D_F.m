% MIRT3D_F Initializes a matrix of cubic B-spline coefficients
% for a single  3D patch. F times the vector of 64 (4x4x4) given control points
% will give a new spatial locations of image voxels within this patch.
%
% Copyright (C) 2007-2010 Andriy Myronenko (myron@csee.ogi.edu)
% also see http://www.bme.ogi.edu/~myron/matlab/MIRT/
% This file is a part of Medical Image Registration Toolbox (MIRT)

function Fi=mirt3D_F(okno)

% create the matrix of weights
  B=[-1 3 -3 1;
    3 -6 3 0;
    -3 0 3 0;
    1 4 1 0]/6;

% create 'okno' points in between of control points
u=linspace(0,1,okno+1)';
u=u(1:end-1);

% create the polynomial matrix 
T=[u.^3 u.^2 u ones(okno,1)];

% compute the precomputed matrix of B-spline basis functions in 1D
B=T*B;

% do kronneker produce to create a 2D matrix
Fi=kron(B,B);

% one more time to create a 3D precomputed matrix
Fi=kron(B,Fi);

