function img_out = tom_xmipp_psd_enhance(img,center,take_log,filter_w1,filter_w2,decay_width,mask_w1,mask_w2)
%TOM_XMIPP_PSD_ENHANCE is a wrapper for xmipp_enhance_psd
%
%   img_out = tom_xmipp_psd_enhance(img,center,take_log,filter_w1,filter_w2,decay_width,mask_w1,mask_w2);
%
%PARAMETERS
%
%  INPUT
%   img                 psd, created with tom_calc_periodogram
%   center              (optional) enable centering (default: true)
%   take_log            (optional) apply log10 to psd (default: true)
%   filter_w1           (optional) Low cut-off frequency for band-pass filtration (default: 0.05)
%   filter_w2           (optional) High cut-off frequency for band-pass filtration (default: 0.2)
%   decay_width         (optional) Decay for the transition bands (default: 0.02)
%   mask_w1             (optional) Inner radius for the annular mask in the frequency domain (default: 0.025)
%   mask_w2             (optional) Outer radius for the annular mask in the frequency domain (default: 0.2)
%
%  OUTPUT
%   img_out             enhanced psd
%
%EXAMPLE
%
%REFERENCES
%
%   http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Enhance_PSD
%
%SEE ALSO
%
%   TOM_CALC_PERIODOGRAM
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
error(nargchk(1, 8, nargin));
error(nargoutchk(1, 1, nargout));

if nargin < 8
    mask_w2 = .2;
end

if nargin < 7
    mask_w1 = .025;
end

if nargin < 6
    decay_width = .02;
end

if nargin < 5
    filter_w2 = .2;
end

if nargin < 4
    filter_w1 = .05;
end

if nargin < 3
    take_log = true;
end

if nargin < 2
    center = true;
end


%cast the input variables to the correct type and call the xmipp function
%img = log(abs(fft2(img)).^2);
img = double(img);

center = logical(center);
take_log = logical(take_log);
filter_w1 = single(filter_w1);
filter_w2 = single(filter_w2);
decay_width = single(decay_width);
mask_w1 = single(mask_w1);
mask_w2 = single(mask_w2);

img_out = tom_xmipp_psd_enhance_wrapper(img,center,take_log,filter_w1,filter_w2,decay_width,mask_w1,mask_w2);