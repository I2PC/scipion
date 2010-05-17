function mask = tom_xmipp_mask(msize,type,origin,varargin)
%TOM_XMIPP_MASK is a wrapper for xmipp_mask
%
%   img_out = tom_xmipp_mask(img,type,varargin);
%
%PARAMETERS
%
%  INPUT
%   msize       size of mask
%   type        type of mask
%   origin      ...
%  
%  OUTPUT
%   data		...
%
%EXAMPLE
%
%REFERENCES
%
%    http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Mask
%
%SEE ALSO
%
%   created by AK 10/11/07
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom
%
%
%check for correct number of input and output arguments
error(nargoutchk(1, 1, nargout));

msizesz = length(msize);
if msizesz > 3 || msizesz < 2
    error('mask size incorrect');
end

if msizesz == 2
    msize(3) = 1;
end

mode = 2;
R1 = 0;
R2 = 0;
pix_width = 0;
H = 0;
sigma = 0;
omega = 0;
rectdim = [0 0 0];
if nargin < 3 || isempty(origin)
    origin = [ceil(msize(1)/2) ceil(msize(2)/2) ceil(msize(3)/2)];
end
s = [0 0];
switch lower(type)
    case 'circular'
        type = 1;
        R1 = varargin{1};
        if R1 < 0
            R1 = abs(R1);
            mode = 1;
        end
    case 'rectangular'
        type = 4;
        rectdim = [varargin{1} varargin{2}];
        if length(varargin) > 2
            rectdim(3) = varargin{3};
        end
        if rectdim(1) < 0 && rectdim(2) < 0 && rectdim(3) < 0
            rectdim = abs(rectdim);
            mode = 1;
        end
    case 'crown'
        type = 2;
        R1 = varargin{1};
        R2 = varargin{2};
        if R1 < 0 && R2 < 0
            R1 = abs(R1);
            R2 = abs(R2);
            mode = 1;
        end
    case 'cylinder'
        type = 3;
        R1 = varargin{1};
        H = varargin{2};
        if R1 < 0 && R2 < 0
            R1 = abs(R1);
            H = abs(H);
            mode = 1;
        end
    case 'cone'
        type = 13;
        R1 = varargin{1};
        if R1 < 0
            R1 = abs(R1);
            mode = 1;
        end
    case 'wedge'
        type = 14;
        R1 = varargin{1};
        R2 = varargin{2};
        if R1 < 0 && R2 < 0
            R1 = abs(R1);
            R2 = abs(R2);
            mode = 1;
        end
    case 'gaussian'
        type = 5;
        sigma = varargin{1};
        if sigma < 0
            sigma = abs(sigma);
            mode = 1;
        end
    case 'raised_cosine'
        type = 6;
        R1 = varargin{1};
        R2 = varargin{2};
        if R1 < 0 && R2 < 0
            R1 = abs(R1);
            R2 = abs(R2);
            mode = 1;
        end
    case 'raised_crown'
        type = 11;
        R1 = varargin{1};
        R2 = varargin{2};
        pix_width = varargin{3};
        if R1 < 0 && R2 < 0
            R1 = abs(R1);
            R2 = abs(R2);
            mode = 1;
        end
    case 'blackman'
        type = 7;
        mode = 1;
    case 'sinc'
        type = 8;
        omega = varargin{1};
        if omega < 0
            omega = abs(omega);
            mode = 1;
        end
        
    otherwise
        error('mask type unknown');
end
        
%cast the input variables to the correct type and call the xmipp function
type = int32(type);
mode = int32(mode);
R1 = double(R1);
R2 = double(R2);
pix_width = double(pix_width);
H = double(H);
sigma = double(sigma);
omega = double(omega);
rectdim = int32(rectdim);
origin = double(origin);
s = int32(s);
msize = double(ones(msize(1),msize(2),msize(3)));
mask = tom_xmipp_mask_wrapper(msize,type,mode,R1,R2,pix_width,H,sigma,omega,rectdim,origin,s);
