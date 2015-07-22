% MIRT3D_grad computes the value of the similarity measure and its gradient

% Copyright (C) 2007-2010 Andriy Myronenko (myron@csee.ogi.edu)
% also see http://www.bme.ogi.edu/~myron/matlab/MIRT/
%
% This file is part of the Medical Image Registration Toolbox (MIRT).

function  [f, Gr, imsmall]=mirt3D_grad(X,  main)

% find the dense transformation for a given position of B-spline control
% points (X).
[Xx,Xy,Xz]=mirt3D_nodes2grid(X, main.F, main.okno);

% Compute the similarity function value (f) and its gradient (dd) at Xx, Xy
% (densely)
[f,ddx,ddy, ddz, imsmall]=mirt3D_similarity(main, Xx, Xy,Xz);

% Find the gradient at each B-spline control point
Gr=mirt3D_grid2nodes(ddx, ddy, ddz, main.F, main.okno, main.mg, main.ng, main.kg);





                        