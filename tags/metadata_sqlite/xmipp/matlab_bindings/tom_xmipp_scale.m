function img_out = tom_xmipp_scale(img,outsize,gridding)
%TOM_XMIPP_SCALE is a wrapper for xmipp_scale
% 
%   img_out = tom_xmipp_scale(img,outsize,gridding);
%
%PARAMETERS
%
%  INPUT
%   img                 2D image or 3D volume
%   outsize             2 or 3 element vector with output dimensions
%   gridding            true or false 
%
%  OUTPUT
%   img_out             scaled image or volume
%
%EXAMPLE
%   img_out = tom_xmipp_scale(img,[64 64],true);
%
%REFERENCES
%
%   http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Scale
%
%SEE ALSO
%   
%
%   created by AK 30/08/07
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
error(nargchk(2, 3, nargin));
error(nargoutchk(1, 1, nargout));

nd = ndims(img);
if  nd ~= 2 && nd ~=3 
    error('only 2D or 3D inputs supported.');
end

if nargin < 3
    gridding = false;
end

osz = length(outsize);
if osz  < 2 || osz > 3
    error('output size argument must be a 2 or 3 element vector');
end

if nd == 2
    outsize(3) = 1;
end

%cast the input variables to the correct type and call the xmipp function
img = double(img);
outsize = double(outsize);
gridding = logical(gridding);
img_out = tom_xmipp_scale_wrapper(img,outsize,gridding);
