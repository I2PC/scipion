function img_out = tom_xmipp_mirror(img,flipstring)
%TOM_XMIPP_MIRROR is a wrapper for xmipp_mirror
%
%   img_out = tom_xmipp_mirror(img,flipstring);
%
%PARAMETERS
%
%  INPUT
%   img                 2D image or 3D volume
%   flipstring          image will be mirrored around these axes: 'x' 'y' or 'xy'
%
%  OUTPUT
%   img_out             mirrored image or volume
%
%EXAMPLE
%   img_out = tom_xmipp_mirror(img,'xy');
%
%REFERENCES
%
%   http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Mirror
%
%SEE ALSO
%   
%
%   created by AK 10/10/07
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
error(nargchk(2, 2, nargin));
error(nargoutchk(1, 1, nargout));

nd = ndims(img);
if  nd ~= 2 && nd ~=3 
    error('only 2D or 3D inputs supported.');
end

flipflag = zeros(3,1,'int8');
if strfind(flipstring,'x')
    flipflag(1) = 1;
end

if strfind(flipstring,'y')
    flipflag(2) = 1;
end

if nd == 3 && strfind(flipstring,'z')
    flipflag(3) = 0;
end

if sum(flipflag,1) == 0
    img_out = img;
    return;
end

%cast the input variables to the correct type and call the xmipp function
img = double(img);
img_out = tom_xmipp_mirror_wrapper(img,flipflag);
%img_out = img_out';