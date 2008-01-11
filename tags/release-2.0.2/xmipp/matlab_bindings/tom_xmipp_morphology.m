function img_out = tom_xmipp_morphology(img,operation,neig,ksize,count)
%TOM_XMIPP_MORPHOLOGY is a wrapper for xmipp_morphology
%
%   img_out = tom_xmipp_morphology(img,operation,neig,ksize,count);
%
%PARAMETERS
%
%  INPUT
%   img                 2D image or 3D volume
%   operation           'dilation','erosion','opening','closing'
%   neig                (optional) for 2D valid neighbourhoods are 4 or 8 (by default), and for 3D: 6, 18 (by default) or 26
%   ksize               (optional) Size of the structuring element (box) (default 1)
%   count               (optional) Minimum number of neighbours with a distinct value to apply the filter
%
%  OUTPUT
%   img_out             morphologically modified image or volume
%
%EXAMPLE
%   img_out = tom_xmipp_morphology(img,'dilation');
%
%REFERENCES
%
%   http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Morphology
%
%SEE ALSO
%   
%
%   created by AK 11/10/07
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
error(nargchk(2, 5, nargin));
error(nargoutchk(1, 1, nargout));

nd = ndims(img);
if  nd ~= 2 && nd ~=3 
    error('only 2D or 3D inputs supported.');
end

if nargin < 5
    count = 0;
end

if nargin < 4
    ksize = 1;
end

if nargin < 3
    if nd == 2
        neig = 8;
    else
        neig = 18;
    end
else
    if (nd == 2 && (neig ~= 4 && neig ~= 8)) || (nd == 3 && (neig ~= 6 && neig ~= 18 && neig ~= 26))
        error('neig parameter is invalid');
    end
end

switch lower(operation)
    case 'dilation'
        op = 1;
    case 'erosion'
        op = 2;
    case 'opening'
        op = 3;
    case 'closing'
        op = 4;
    otherwise
        error('Operation not supported.');
end

        
%cast the input variables to the correct type and call the xmipp function
img = double(img);
op = int32(op);
ksize = int32(ksize);
count = int32(count);
neig = int32(neig);
img_out = tom_xmipp_morphology_wrapper(img,op,ksize,count,neig);