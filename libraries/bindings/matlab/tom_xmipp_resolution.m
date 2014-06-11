function st = tom_xmipp_resolution(img,ref,objectpixelsize)
%TOM_XMIPP_RESOLUTION is a wrapper for xmipp_resolution
%
%   st = tom_xmipp_resolution(img,ref,objectpixelsize);
%
%PARAMETERS
%
%  INPUT
%   img                 2D or 3D image or volume
%   ref                 2D or 3D reference
%   objectpixelsize     object pixel size of the 2 volumes in Angstroems
%  
%  OUTPUT
%   st.freq             frequency
%   st.dpr              differential phase residual
%   st.frc              fourier shell coefficient
%   st.frc_noise
%
%EXAMPLE
%
%
%REFERENCES
%
% http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Resolution
%
%SEE ALSO
%   
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
error(nargchk(3, 3, nargin))
error(nargoutchk(1, 1, nargout));

sz_img = size(img);
sz_ref = size(ref);

if sum(sz_img-sz_ref) > 0
    error('image and reference must be of the same size');
end

if objectpixelsize <= 0
    error('object pixel size must be positive.');
end

        
%cast the input variables to the correct type and call the xmipp function
img = double(img);
ref = double(ref);
objectpixelsize = single(objectpixelsize);
st = tom_xmipp_resolution_wrapper(img,ref,objectpixelsize);