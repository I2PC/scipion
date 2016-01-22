% MIRT3D_REGISTER The main function for non-rigid registration 
% of a single pair of 3D images using cubic B-spline
% based transformation parametrization. 
%
%   Input
%   ------------------ 
%   refim       Reference (or fixed) 2D image with intensities [0..1]
%   im          Source (or float) 2D image to be deformed with intensities [0..1]
%               If you want to exlude certain image parts from the
%               computation of the similarity measure, set those parts as
%               NaNs (masking).
%
%   main        a structure of main MIRT options
%
%       .similarity=['SSD','CC','SAD','RC', 'CD2', 'MS','MI'] (default SSD) 
%                    SSD - Sum of squared differences 
%                    CC  - Correlation Coefficient 
%                    SAD - Sum of absolute differences
%                    RC -  Residual Complexity (monomodal, nonstationary slow-varying intensity distortions)
%                    CD2, MS - Ultrasound similarity measures
%                    MI - Mutual Inormation (multimodal)
%               Similarity measure. All current similarity measures are
%               listed in mirt3D_similarity. You can easily add your own
%               one there. If you do please email me and I'll include it in MIRT
%               with full acknowledgment.
%       .okno (default 16) mesh window size between the B-spline
%               control points. The mesh cell is square. The smaller the
%               window the more complex deformations are possible, but
%               also more regularization (main.lambda) is required. 
%       .subdivide (default 3) - a number of hierarchical levels. 
%               E.g. for main.subdivide=3 the registration is carried sequentially
%               at image size 2^(3-1)=4 times smaller, 2^(2-1)=2 times smaller and the 
%               original size. The mesh window size remain the same for all levels.
%       .lambda (default 0.01) - a regularization weight. A regularization
%               is defined as Laplacian (or curvature) penalization of the
%               displacements of B-spline control points. Set main.single=1 to see the mesh,
%               if your mesh is getting too much deformed or even unnaturally folded,
%               set lambda to higher values.(try [0..0.4]). See my thesis 5.7
%       .alpha  (default 0.1) - a parameter of the similarity measure, for
%               e.g. alpha value of the Residual Complexity (try [0.01..1]) or scaling (parameter D) in CD2 and MS.
%       .single (default 0) show mesh deformation at each iteration
%               0 - do not show, 1 or anything else - show.
%       .ro ([0..0.99]) a correlation parameter of MS similarity measure
%       .MIbins (default 64)  Number of bins to use for the MI similarity measure
%
%   optim       a structure of optimization options
%
%       .maxsteps (default 300) Maximum number of iterations. If at the
%               final iterations, the similarity measure is still getting
%               significantly decreased, you may need to set more maximum
%               iterations.
%       .fundif (default 1e-6) Tolerance, stopping criterion. 
%       .gamma (default 1) Important parameter: Initial optimization step size.
%               During optimization the optimizer is taking
%               steps starting from (optim.gamma) to decrease
%               the similarity measure. If the step is too big
%               it will be adjusted to smaller as gamma=gamma*anneal (see
%               below). During the registration look at the
%               value of gamma, if it stays constant at an
%               initial level for too long, it means your
%               initial step is too small, and the optimization
%               will take longer. Ideally gamma is constant
%               during first several iterations and then slowly
%               decreasing.
%       .anneal (default 0.9) The multiplicative constant to update the
%               step size
%
%   Output
%   ------------------ 
%   res         a structure of the resulting parameters:
%
%       .X      a 4D matrix of final B-spline control points positions
%               The 4th dimension size is 3, it indicates the coordinate
%               x,y and z.
%       .okno   window size/spacing between the control points (equal to
%               the initial main.okno).
%   im_int      Deformed float image (result)
%
%           
%
%   Examples
%   --------
%
%   See many detailed examples in the 'examples' folder.
%
%   See also mirt3D_registration, mirt2D_register

% Copyright (C) 2007-2010 Andriy Myronenko (myron@csee.ogi.edu)
% also see http://www.bme.ogi.edu/~myron/matlab/MIRT/
%
%     This file is part of the Medical Image Registration Toolbox (MIRT).
%
%     The source code is provided under the terms of the GNU General Public License as published by
%     the Free Software Foundation version 2 of the License.

function [res, im_int]=mirt3D_register(refim, im, main, optim)

% Check the input options and set the defaults
if nargin<2, error('mirt3D_register error! Not enough input parameters.'); end;
% Check the proper parameter initilaization
[main,optim]=mirt_check(main,optim,nargin);

% checking for the possible errors
if numel(size(im))~=numel(size(refim)), error('The dimensions of images are not the same.'); end;
if size(im,1)~=size(refim,1), error('The images must be of the same size.'); end;
if size(im,2)~=size(refim,2), error('The images must be of the same size.'); end;
if size(im,3)~=size(refim,3), error('The images must be of the same size.'); end;


tic
disp('MIRT: Starting 3D registration...');
disp('Using the following parameters:');
main
optim

% Original image size 
dimen=size(refim);
% Size at the smallest hierarchical level, when we resize to smaller
M=ceil(dimen/2^(main.subdivide-1));
% Generate B-splines mesh of control points equally spaced with main.okno spacing
% at the smallest hierarchical level (the smallest image size)
[x, y, z]=meshgrid(1-main.okno:main.okno:M(2)+2*main.okno, 1-main.okno:main.okno:M(1)+2*main.okno, 1-main.okno:main.okno:M(3)+2*main.okno);

% the size of the B-spline mesh is in (main.mg, main.ng, main.kg) at the smallest
% hierarchival level.
[main.mg, main.ng, main.kg]=size(x); 

% new image size at the smallest hierarchical level
% this image size is equal or bigger than the original resized image at the
% smallest level. This size includes the image and the border of nans when the image size can not be
% exactly divided by integer number of control points.
main.siz=[(main.mg-3)*main.okno (main.ng-3)*main.okno (main.kg-3)*main.okno];

main.X=cat(4,x,y,z);  % Put x, y and z control point positions into a single mg x ng x kg x 3  4Dmatrix
main.Xgrid=main.X;    % save the regular grid (used for regularization)

main.F=mirt3D_F(main.okno); % Init B-spline coefficients

% Leave only the image described by control points.
% At the original (large) image size calculate
% the size described my B-spline control points,
% which will be larger or equal to the original image size
% and initialize with NaNs. Then patch with the original images,
% so that the original images (refim, im) now include the border of NaNs.
% NaNs here signalizes the values to be ignored during the registration
tmp=nan(2^(main.subdivide-1)*main.siz);
tmp(1:dimen(1),1:dimen(2),1:dimen(3))=refim;  refim=tmp;
tmp(1:dimen(1),1:dimen(2),1:dimen(3))=im;     im=tmp;
clear tmp;

% Go across sublevels
for level=1:main.subdivide
    
    % update the size of B-spline mesh to twice bigger 
    % only do it for the 2nd or higher levels 
    if level>1,
        main.mg=2*main.mg-3; % compute the new mesh size
        main.ng=2*main.ng-3;
        main.kg=2*main.kg-3;
    end
    main.level=level;
    
    % current image size
    main.siz=[(main.mg-3)*main.okno (main.ng-3)*main.okno (main.kg-3)*main.okno];
    main.K=mirt3D_initK([main.mg main.ng main.kg]); % Init Laplacian eigenvalues (used for regularization)
    
    main.refimsmall=mirt3D_imresize(refim,main.siz); % resize images
    main.refimsmall(main.refimsmall<0)=0;  
    main.refimsmall(main.refimsmall>1)=1;
    
    imsmall=mirt3D_imresize(im,main.siz);
    imsmall(imsmall<0)=0;
    imsmall(imsmall>1)=1;
    
    [gradx, grady, gradz]=gradient(imsmall);
    main.imsmall=cat(4,imsmall, gradx,grady, gradz); % concatenate the image with its gradients
    % main.imsmall is a 4D matrix - a concatination of the image and its x,
    % y and z gradients at the given hierarchical level. The reason to
    % concatinate them together is to save time later on image interpolation
    % using mirt3D_mexinterp. The image and its gradients have to be
    % interpolated at the same coordinates, which can be done more
    % efficiently then interpolating each of them individially.
    
    
    % a single level 3D non-rigid image registration
    [main.X, result]=mirt3D_registration(main.X, main, optim);
    
    % if the sublevel is not last prepare for the next level
    if level<main.subdivide,
        main.X=mirt3D_subdivide(main.X, 1);
        main.Xgrid=mirt3D_subdivide(main.Xgrid,  1);
    end;
    
end

% Prepare the output
res.X=main.X;
res.Xgrid=main.Xgrid;
res.okno=main.okno;

% because we have appended the initial images with the border of NaNs during the
% registration, now we want to remove that border and get the result of the
% initial image size
im_int=zeros(dimen); [M,N,K]=size(result);
im_int(1:min(dimen(1),M),1:min(dimen(2),N),1:min(dimen(3),K))=result(1:min(dimen(1),M),1:min(dimen(2),N),1:min(dimen(3),K));


disp('MIRT: 3D non-rigid registration is successfully completed.')
disptime(toc);



