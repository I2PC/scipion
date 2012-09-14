function img_out = tom_xmipp_scale_pyramid(img,operation,levels)
%TOM_XMIPP_SCALE_PYRAMID is a wrapper for xmipp_scale_pyramid
%
%   img_out = tom_xmipp_scale_pyramid(img,operation,levels);
%
%PARAMETERS
%
%  INPUT
%   img                 2D image or 3D volume
%   operation           'expand' or 'reduce'
%
%  OUTPUT
%   img_out             scaled image
%
%EXAMPLE
%   img_out = tom_xmipp_scale_pyramid(img,'reduce',2);
%
%REFERENCES
%
%   http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Scale_pyramid
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

if nargin < 2
    levels = 1;
end

if strcmpi(operation,'expand') == true
    operation = true;
else
    operation = false;
end
        

%cast the input variables to the correct type and call the xmipp function
img = double(img);
levels = int32(levels);
img_out = tom_xmipp_scale_pyramid_wrapper(img,operation,levels);