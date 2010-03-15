function img_out = tom_xmipp_scale(img,angs,axis,align_z,gridding,wrap)
%TOM_XMIPP_ROTATE is a wrapper for xmipp_rotate
%
%   img_out = tom_xmipp_rotate(img,angs,axis,align_z,gridding,wrap);
%
%PARAMETERS
%
%  INPUT
%   img                 2D image or 3D volume
%   angs                rotation angle for 2D, [rot tilt psi] for 3D
%   axis                (optional) 3 element vector rotation axis or empty
%   align_Z             (optional) 3 element vector of point to align to Z axis or empty
%   gridding            (optional) true or false (default: false)
%   wrap                (optional) true or false (default: true)
%
%  OUTPUT
%   img_out             rotated image or volume
%
%EXAMPLE
%   img_out = tom_xmipp_rotate(img,30);
%
%REFERENCES
%
%   http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Rotate
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
error(nargchk(2, 6, nargin));
error(nargoutchk(1, 1, nargout));

nd = ndims(img);

if  nd ~= 2 && nd ~=3 
    error('only 2D or 3D inputs supported.');
end

if nargin < 6
    wrap = true;
end
if nargin < 5
    gridding = false;
end

if nargin < 4
    align_z = [];
end

if nargin < 3 || isempty(axis)
    axis = [0 0 1];
end

if length(axis) ~= 3
    error('axis argument must be a 3 element vector');
end

tform = zeros(3,3);
angsz = size(angs);
if angsz(2) == 3 && angsz(1) == 1
    mode = 1;
elseif angsz(2) == 3 && angsz(1) == 3
    mode = 4;
    tform = angs;
else
    if ~isempty(align_z);
        mode = 2;
        axis = align_z;
        angs = [0 0 0];
    else
        mode = 3;
    end
end

if nd == 2
    angs(2:3) = 0;
end


%cast the input variables to the correct type and call the xmipp function
img = double(img);
angs = double(angs);
axis = double(axis);
gridding = logical(gridding);
wrap = logical(wrap);
tform = double(tform);
img_out = tom_xmipp_rotate_wrapper(img,angs,axis',mode,gridding,wrap,tform);