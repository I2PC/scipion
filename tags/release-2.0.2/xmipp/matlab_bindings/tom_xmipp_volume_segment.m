function [seg_mask, vol_seg] = tom_xmipp_volume_segment(vol,sampling,mass,type,enable_threshold,threshold,wang_radius,probabilistic)
%TOM_XMIPP_VOLUME_SEGMENT is a wrapper for xmipp_volume_segment
%
%   vol_seg = tom_xmipp_volume_segment(vol,sampling,mass,type,enable_threshold,threshold,wang_radius,probabilistic);
%
%PARAMETERS
%
%  INPUT
%   vol                 3D volume
%   sampling            sampling rate [Angstroems/voxel]
%   mass                value to segment
%   type                'voxels','daltons','amino acids'
%   enable_threshold    (optional) enable threshold value (default: false)
%   threshold           (optional) threshold value, ignored if enable_threshold is false
%   wang_radius         (optional) Radius [pix] for B.C. Wang cone  (default: 3)
%   probabilistic       (optional) Calculate probabilistic solvent mask (default: false)
%
%  OUTPUT
%   seg_mask             segmentation mask
%   vol_seg              (optional) segmented volume
%
%EXAMPLE
%   
%
%REFERENCES
%
%   http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Segment
%
%SEE ALSO
%   
%
%   created by AK 1/11/07
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
error(nargchk(4, 8, nargin));
error(nargoutchk(1, 2, nargout));

nd = ndims(vol);
if   nd ~=3 
    error('only 3D inputs supported.');
end

if nargin < 8
    probabilistic = 0;
end

if nargin < 7
    wang_radius = 3;
end

if nargin < 6
    threshold = 0;
end

if nargin < 5
    enable_threshold = 0;
end

switch type
    case 'voxels'
        type = 1;
    case 'daltons'
        type = 2;
    case 'amino acids'
        type = 3;
    otherwise
        error('segmentation type not supported.');
end


%cast the input variables to the correct type and call the xmipp function
vol = double(vol);
sampling = double(sampling);
mass = double(mass);
type = int32(type);
enable_threshold = int32(enable_threshold);
threshold = double(threshold);
wang_radius = int32(wang_radius);
probabilistic = int32(probabilistic);
seg_mask = tom_xmipp_volume_segment_wrapper(vol,sampling,mass,type,enable_threshold,threshold,wang_radius,probabilistic);
if nargout == 2
    vol_seg = seg_mask .* vol;
end