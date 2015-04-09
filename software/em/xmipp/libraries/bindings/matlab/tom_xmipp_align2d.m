function st = tom_xmipp_align2d(img,ref,mode,max_shift,max_rot,psi_interval,Rin,Rout,outside)
%TOM_XMIPP_ALIGN2D is a wrapper for xmipp_align2d
%
%   img_out = tom_xmipp_align2d(img,ref,mode,max_shift,max_rot,psi_interval,Rin,Rout,outside);
%
%PARAMETERS
%
%  INPUT
%   img                 2D image
%   ref                 2D reference
%   mode                (optional) 'rot' 'trans' 'complete' (default: complete)
%   max_shift           (optional) maximum shift (default 0)
%   max_rot             (optional) maximum rotation (default 0)
%   psi_interval        (optional) angular increment for complete search mode
%   Rin                 (optional) only the region between radii Ri and Ro (in pixels) will be considered for rotational correlation
%   Rout
%   outside             (optional) values outside the mask will be replaced by this value (default: 0)
%
%  OUTPUT
%
%EXAMPLE
%
%   st = tom_xmipp_align2d(img,ref);
%
%REFERENCES
%
%   http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Align2D
%
%SEE ALSO
%   
%
%   created by AK 31/08/07
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
error(nargchk(1, 9, nargin));
error(nargoutchk(1, 1, nargout));

if nargin < 9
    outside = 0;
end

if nargin < 8
    Rin = 0;
end

if nargin < 7
    Rout = size(img,1)./2;
end

if nargin < 6
    psi_interval = 10;
end

if nargin < 5
    max_rot = 0;
end

if nargin < 4
    max_shift = 0;
end

if nargin < 3
	mode = 'complete';
end

switch mode 
    case 'trans'
        mode = 1;
    case 'rot'
        mode = 2;
    case 'complete'
        mode = 3;
end

%cast the input variables to the correct type and call the xmipp function
img = double(img);
ref = double(ref);
mode = int32(mode);
max_shift = single(max_shift);
max_rot = single(max_rot);
psi_interval = single(psi_interval);
Rin = single(Rin);
Rout = single(Rout);
outside = double(outside);

st = tom_xmipp_align2d_wrapper(img,ref,mode,max_shift,max_rot,psi_interval,Rin,Rout,outside);