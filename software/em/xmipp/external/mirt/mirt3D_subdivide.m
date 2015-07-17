% This function subdivide the current mesh of B-spline control 
% points to one level up (twise denser). The new size will be
% mg=2*mg-3. The function can subdivide a single mesh or a sequence
% meshes simultaneusly.
%
% The idea is to take each 3D cube of 8 (2x2x2) of B-spline control points
% and upsample it to form a new 3D cube of 27 (3x3x3). And repeat for
% for all B-spline control points
% After upsampling is done, we also need to remove the control points
% around the image border, because they are useless at the higher
% resolution level, and will only add to the computational time.

% Copyright (C) 2007-2010 Andriy Myronenko (myron@csee.ogi.edu)
% also see http://www.bme.ogi.edu/~myron/matlab/MIRT/
%
% This file is part of the Medical Image Registration Toolbox (MIRT).

function Y=mirt3D_subdivide(X, M)

[mg,ng, kg, tmp, tmp]=size(X);

x=squeeze(X(:,:,:,1,:));
y=squeeze(X(:,:,:,2,:));
z=squeeze(X(:,:,:,3,:));

xnew=zeros(mg, 2*ng-2, kg, M);
ynew=zeros(mg, 2*ng-2, kg, M);
znew=zeros(mg, 2*ng-2, kg, M);

xfill=(x(:,1:end-1,:,:)+x(:,2:end,:,:))/2;
yfill=(y(:,1:end-1,:,:)+y(:,2:end,:,:))/2;
zfill=(z(:,1:end-1,:,:)+z(:,2:end,:,:))/2;


for i=1:ng-1
   xnew(:,2*i-1:2*i,:,:)=cat(2,x(:,i,:,:), xfill(:,i,:,:)); 
   ynew(:,2*i-1:2*i,:,:)=cat(2,y(:,i,:,:), yfill(:,i,:,:)); 
   znew(:,2*i-1:2*i,:,:)=cat(2,z(:,i,:,:), zfill(:,i,:,:));
end

x=xnew(:,2:end,:,:);
y=ynew(:,2:end,:,:); 
z=znew(:,2:end,:,:); 
%

xnew=zeros(2*mg-2, 2*ng-3, kg, M);
ynew=zeros(2*mg-2, 2*ng-3, kg, M);
znew=zeros(2*mg-2, 2*ng-3, kg, M);


xfill=(x(1:end-1,:,:,:)+x(2:end,:,:,:))/2;
yfill=(y(1:end-1,:,:,:)+y(2:end,:,:,:))/2;
zfill=(z(1:end-1,:,:,:)+z(2:end,:,:,:))/2;

for i=1:mg-1
   ynew(2*i-1:2*i,:,:,:)=cat(1,y(i,:,:,:), yfill(i,:,:,:)); 
   xnew(2*i-1:2*i,:,:,:)=cat(1,x(i,:,:,:), xfill(i,:,:,:)); 
   znew(2*i-1:2*i,:,:,:)=cat(1,z(i,:,:,:), zfill(i,:,:,:));
end

x=xnew(2:end,:,:,:);
y=ynew(2:end,:,:,:);
z=znew(2:end,:,:,:);
%
xnew=zeros(2*mg-3, 2*ng-3, 2*kg-2, M);
ynew=zeros(2*mg-3, 2*ng-3, 2*kg-2, M);
znew=zeros(2*mg-3, 2*ng-3, 2*kg-2, M);


xfill=(x(:,:,1:end-1,:)+x(:,:,2:end,:))/2;
yfill=(y(:,:,1:end-1,:)+y(:,:,2:end,:))/2;
zfill=(z(:,:,1:end-1,:)+z(:,:,2:end,:))/2;

for i=1:kg-1
   ynew(:,:,2*i-1:2*i,:)=cat(3,y(:,:,i,:), yfill(:,:,i,:)); 
   xnew(:,:,2*i-1:2*i,:)=cat(3,x(:,:,i,:), xfill(:,:,i,:)); 
   znew(:,:,2*i-1:2*i,:)=cat(3,z(:,:,i,:), zfill(:,:,i,:));
end

x=xnew(:,:,2:end,:);
y=ynew(:,:,2:end,:);
z=znew(:,:,2:end,:);

Y=cat(5, x, y, z);
Y=permute(Y,[1 2 3 5 4]);
Y=2*Y-1;