function st = tom_xmipp_adjust_ctf(psd,Dz,voltage,objectPixelSize,ctfmodelSize,Cs,min_freq,max_freq,Ca,enhance_filter_min,enhance_filter_max,enhance_weight)
%TOM_XMIPP_ADJUST_CTF is a wrapper for xmipp_adjust_ctf
%
%   st = tom_xmipp_adjust_ctf(psd,Dz,voltage,objectPixelSize,ctfmodelSize,Cs,min_freq,max_freq,Ca,enhance_filter_min,enhance_filter_max,enhance_weight);
%
%PARAMETERS
%
%  INPUT
%   psd                 periodogram, generated with tom_calc_periodogram
%   Dz                  initial defocus value
%   voltage             accelerating voltage     
%   objectPixelSize     object pixel size
%   ctfmodelSize        (optional) size of the output CTF models (default 0)
%   Cs                  (optional) Cs value (default 2)
%   min_freq            (optional) minimal frequency to consider in the fit (default 0.03)
%   max_freq            (optional) maximal frequency to consider in the fit (default 0.3)
%   Ca                  (optional) Ca value (default 2) 
%   enhance_filter_min  (optional) minimal frequency cutoff of Fourier filter for PSD enhancement (0.02)
%   enhance_filter_max  (optional) maximal frequency cutoff of Fourier filter for PSD enhancement (0.15)
%   enhance_weight      (optional) enhance_weight (default 5)
%
%  OUTPUT
%   st.DeltafU          defocus Value 1
%   st.DeltafV          defocus Value 2
%   st.AzimuthalAngle   angle of the astigmatism
%   st.CTFmodelhalf     half circle CTF model (only if ctfmodelSize > 0)
%   st.CTFmodelquadrant quadrant circle CTF model (only if ctfmodelSize > 0)
%
%
%EXAMPLE
%
%REFERENCES
%   
%   http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Adjust_CTF
%
%SEE ALSO
%   
%   TOM_CALC_PERIODOGRAM, TOM_XMIPP_ENHANCE_PSD,
%   TOM_XMIPP_CTF_CORRECT_PHASE
%
%   created by AK & COS 27/08/07
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
error(nargchk(4, 12, nargin));
error(nargoutchk(1, 1, nargout));

if nargin < 12
    enhance_weight = 5;
end

if nargin < 11
    enhance_filter_max = 0.15;
end

if nargin < 10
    enhance_filter_min = 0.02;
end

if nargin < 9
    Ca = 2;
end

if nargin < 8
    max_freq = 0.3;
end

if nargin < 7
    min_freq = 0.03;
end

if nargin < 6
    Cs = 2;
end

if nargin < 5
    ctfmodelSize = 0;
end

    
%cast the input variables to the correct type and call the xmipp function
psd = double(psd);
Dz = double(Dz);
voltage = double(voltage);
objectPixelSize = double(objectPixelSize);
Cs = double(Cs);
ctfmodelSize = double(ctfmodelSize);
min_freq = double(min_freq);
max_freq = double(max_freq);
Ca = double(Ca);
enhance_filter_min = double(enhance_filter_min);
enhance_filter_max = double(enhance_filter_max);
enhance_weight = double(enhance_weight);

st = tom_xmipp_adjust_ctf_wrapper(psd,min_freq,max_freq,ctfmodelSize,Ca,objectPixelSize,voltage,Cs,Dz,enhance_filter_min,enhance_filter_max,enhance_weight);
