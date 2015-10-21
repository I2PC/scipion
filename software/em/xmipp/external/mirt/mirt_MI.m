% MIRT_MI  computes the negative Mutual Information and its
% voxel-wise gradient.
%
% Input 
% refim - reference image 2D or 3D
% im - template/source/float image, same size as refim
% bins - the number of intensity levels to use
%
% Output
% MI - the value of negative MI
% GMI - an array with voxel-wise (or pixel-wise in 2D case) gradients of
% the similarity measure.

% Copyright (C) 2007-2010 Andriy Myronenko (myron@csee.ogi.edu)
% also see http://www.bme.ogi.edu/~myron/matlab/MIRT/
%
% This file is part of the Medical Image Registration Toolbox (MIRT).

function [MI, GMI]=mirt_MI(refim,im, bins)

% remove nans
refim(isnan(refim))=0;
im(isnan(im))=0;
M=numel(im);

%%Scott's rule for Gaussian width 
%d=2; sigma=(M*(d+2))^(-1/(d+2))*std(refim(:),0);
sigma=0.05;  % default Gaussian width (images must be within (0..1))

% build the segmentation matrices (each column indicates the grayscale level segment)
PI = sparse((1:M)',round(im(:)*(bins-1))+1,true(M,1),M,bins,M);
PJ = sparse((1:M)',round(refim(:)*(bins-1))+1,true(M,1),M,bins,M);

% construct the smoothing matrix
vbins=(0:bins-1)/(bins-1);
G=exp(-bsxfun(@minus,vbins',vbins).^2/(2*sigma^2));
G=G/sqrt(M*(2*pi*sigma^2));

% compute joint and marginal distributions

% If you get an Error here uncomment the line belowe instead
% (For some reason Matlab 7 stopped supporting product of 2 logical matrices)
Pij=G*(PI'*PJ)*G;
%Pij=G*(PI'*double(PJ))*G;

Pi=(sum(PI)*G)'/sqrt(M);
Pj=(sum(PJ)*G)'/sqrt(M);

% negative Mutual Information (-MI)
PiPj=Pi(:)*Pj(:)';
L=log((Pij+(Pij==0))./(PiPj+(PiPj==0)));
MI=-sum(Pij(:).*L(:));

% gradient of MI
%%% GMI=im(:).*full(sum((PI*(G*(L'+1)*G)).*PJ,2))-sum((PI*(G*diag(vbins)*(L'+1)*G)).*PJ,2);
% same as above, but faster
GL=G*(L+1)';
A=(GL*G)*PI';
B=(GL*diag(vbins)*G)*PI';
GMI=(im(:).*A(PJ')-B(PJ'))/sigma^2;

% reshape to the image size
GMI=reshape(GMI,size(im));