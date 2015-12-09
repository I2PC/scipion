% im=mirt3D_imresize(im, newsiz)  resize the 3D image im
% to the new (but smaller) size newsiz. For examples two resize the image to 3 times
% smaller use im=mirt3D_imresize(im,size(im)/3);

% Copyright (C) 2007-2010 Andriy Myronenko (myron@csee.ogi.edu)
% also see http://www.bme.ogi.edu/~myron/matlab/MIRT/
%
%  This file is part of the Medical Image Registration Toolbox (MIRT).


function im=mirt3D_imresize(im, newsiz)


siz=size(im);    % current size

if sum(siz==newsiz)~=3 % check if the new size is not the same

coef=isnan(im);  % if we've used masking through nans, save them
im(coef)=0;      % set nans to zeros

M=siz./newsiz;   % the fraction between current and old sizes 
 

% if resizing to smaller, remove the higher frequency content (low pass filtering through fft)
if sum(M<1)==0,    
    fim=fftshift(fftn(im));
    fim2=zeros(siz);
    range1=round(0.5*(1-1./M).*(siz-1)+1);
    range2=round(0.5*(1+1./M).*(siz-1)+1);
    fim2(range1(1):range2(1),range1(2):range2(2),range1(3):range2(3))=fim(range1(1):range2(1),range1(2):range2(2),range1(3):range2(3));
    im=abs(ifftn(fim2));

    % put nans back
    im(coef)=nan;
end

% prepare the coordinates of the smaller image
x=linspace(1,siz(2),newsiz(2));
y=linspace(1,siz(1),newsiz(1));
z=linspace(1,siz(3),newsiz(3));
[x,y,z]=meshgrid(x,y,z);

% linear interpolate
im=mirt3D_mexinterp(im,x,y,z);

end