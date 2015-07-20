% mirt3D_similarity  The function computes the current similarity measure
% value and its dense (voxel-wise) gradients

% Copyright (C) 2007-2010 Andriy Myronenko (myron@csee.ogi.edu)
% also see http://www.bme.ogi.edu/~myron/matlab/MIRT/
%
% This file is part of the Medical Image Registration Toolbox (MIRT).
function [f,ddx,ddy, ddz,imsmall]=mirt3D_similarity(main, Xx, Xy, Xz)

% interpolate image and its gradients simultaneously
% main.imsmall is a 4D array where
% main.imsmall(:,:,:,1) - is the float image to be deformed
% main.imsmall(:,:,:,2:4) - its x,y,z gradient
im_int=mirt3D_mexinterp(main.imsmall, Xx, Xy, Xz);
clear Xx Xy Xz tmp;

% isolate the interpolated 3D image
imsmall=im_int(:,:,:,1); 

% Compute the similarity function value (f) and its gradient (dd)
switch lower(main.similarity)
   
   % sum of squared differences
   case 'ssd' 
        dd=imsmall-main.refimsmall;
        f=nansum(dd(:).^2)/2;
        
   % sum of absolute differences       
   case 'sad' 
        dd=imsmall-main.refimsmall;
        f=nansum(sqrt(dd(:).^2+1e-10));
        dd=dd./sqrt(dd.^2+1e-10);     
    
   % correlation coefficient     
   case 'cc' 
       
       SJ=main.refimsmall-nansum(main.refimsmall(:))/numel(main.refimsmall);
       SI=imsmall-nansum(imsmall(:))/numel(imsmall);
       
       a = nansum(imsmall(:).*SJ(:))/nansum(imsmall(:).*SI(:));
       f=-a*nansum(imsmall(:).*SJ(:));
       dd=-2*(a*SJ-a^2*SI);     
        
   % Residual Complexity: A. Myronenko, X. Song: "Image Registration by
   % Minimization of Residual Complexity.", CVPR'09
   case 'rc' 
      
        rbig=imsmall-main.refimsmall;
        
        [y,x]=find_imagebox(rbig); r=rbig(y,x);
        r(isnan(r))=nanmean(r(:));

        Qr=mirt_dctn(r);
        Li=Qr.^2+main.alpha;

        f=0.5*sum(log(Li(:)/main.alpha));
        
        r=mirt_idctn(Qr./Li);
        dd=zeros(size(rbig)); 
        dd(y,x)=r;
     
   % CD2 similarity measure: Cohen, B., Dinstein, I.: New maximum likelihood motion estimation schemes for
   % noisy ultrasound images. Pattern Recognition 35(2),2002
   case 'cd2' 
        
        f=(imsmall-main.refimsmall)/main.alpha;
        dd=2*tanh(f);
        f=2*nansum(log(cosh(f(:))));
        
   % MS similarity measure: Myronenko A., Song X., Sahn, D. J. "Maximum Likelihood Motion Estimation
   % in 3D Echocardiography through Non-rigid Registration in Spherical Coordinates.", FIMH 2009
   case 'ms' 
        f=(imsmall-main.refimsmall)/main.alpha;
        coshd2=cosh(f).^2;
        dd=tanh(f).*(2*coshd2+main.ro)./(coshd2-main.ro);
        f=nansum(1.5*log(coshd2(:)-main.ro)-0.5*log(coshd2(:)));
        
  % (minus) Mutual Information: Paul A. Viola "Alignment by Maximization of Mutual Information"   
  case 'mi' 
        % MI computation is somewhat more involved, so let's compute it in separate function
        [f, dd]=mirt_MI(main.refimsmall,imsmall,main.MIbins);

    otherwise
        error('Similarity measure is wrong. Supported values are RC, SAD, SSD, CD2, MS')
end


% Multiply by interpolated image gradients
ddx=dd.*im_int(:,:,:,2); ddx(isnan(ddx))=0;
ddy=dd.*im_int(:,:,:,3); ddy(isnan(ddy))=0;
ddz=dd.*im_int(:,:,:,4); ddz(isnan(ddz))=0;

% This subfunctions finds the coordinates of the largest square
% within the image that has no NaNs (not affected by interpolation)
% It can be useful to ignore the border artifacts caused by interpolation,
% or when the actual image has some black border around it, that you don't
% want to take into account.
function [y,x,z]=find_imagebox(im)

ind=find(~isnan(im));
[i,j,k]=ind2sub(size(im),ind);

n=0; % extra border size to cut
y=min(i)+n:max(i)-n;
x=min(j)+n:max(j)-n;
z=min(k)+n:max(k)-n;
