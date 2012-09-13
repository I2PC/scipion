function img_out = tom_xmipp_normalize(img,method,mask)
%TOM_XMIPP_NORMALIZE is a wrapper for xmipp_normalize
%
%   img_out = tom_xmipp_normalize(img,method,mask);
%
%PARAMETERS
%
%  INPUT
%   img                 2D image
%   method              normalization method, valid choices are: 
%                       'OLDXMIPP','NEAR_OLDXMIPP','NEWXMIPP','MICHAEL','NEWXMIPP2','RAMP'
%   mask                (optional) mask to use with normalization  
%
%  OUTPUT
%   img_out             normalized image
%
%EXAMPLE
%   mask = tom_sphere(size(img),size(img,1)./2,2);
%   img_out = tom_xmipp_normalize(img,'Ramp',mask);
%
%REFERENCES
%
%   http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Normalize
%
%SEE ALSO
%   
%
%   created by AK & COS 24/08/07
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

%check for correct number of input and output arguments
error(nargchk(2, 4, nargin));
error(nargoutchk(1, 1, nargout));

if nargin < 3
    mask = ones(size(img),'int32');
end

%check if image and mask are of the same dimensions
if sum(size(img)==size(mask)) < 2
    error('Image and mask are of unequal size');
end

%check if image and mask are 2D
if ndims(img) ~= 2 || ndims(mask) ~= 2
    error('Image and mask must be 2D arrays');
end

%convert the method into an integer
switch lower(method)
%    case 'none'
%        method = 0;
    case 'oldxmipp'
        method = 1;
    case 'newxmipp'
        method = 3;
    case 'michael'
        method = 4;
    case 'newxmipp2'
        method = 5;
%    case 'random'
%        method = 6;
    case 'ramp' 
        method = 7;
%    case 'neighbour'
%        method = 8;
    otherwise
        error('Not a valid xmipp normalization method');
end

%cast the input variables to the correct type and call the xmipp function
img = double(img);
mask = int32(mask);
method = int32(method);
img_out = tom_xmipp_normalize_wrapper(img,mask,method);