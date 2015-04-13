function img_out = tom_xmipp_ctf_correct_phase(img,st,method,epsilon)
%TOM_XMIPP_CTF_CORRECT_PHASE is a wrapper for xmipp_ctf_correct_phase
%
%   img_out = tom_xmipp_ctf_correct_phase(img,st,method,epsilon);
%
%PARAMETERS
%
%  INPUT
%   img                 2D image 
%   st                  CTF model structure (output of tom_xmipp_adjust_ctf)
%   method (optional)   'remove','leave','divide' (default: leave)
%                         * remove  set all small values to 0. Correct phase of big values
%                         * leave leave all small values as they are. Correct phase of big values
%                         * divide correct amplitude and phase by dividing by the CTF. 
%                           If the value is small then leave it as it is. 
%   epsilon (optional)   CTF values whose absolute value is smaller than
%                        epsilon are assumed to be small (default: 0)
%
%  OUTPUT
%   img_out             phase corrected image
%
%EXAMPLE
%
%REFERENCES
%
%   http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Correct_Phase
%
%SEE ALSO
%
%   TOM_XMIPP_ADJUST_CTF
%
%   created by AK 22/10/07
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
error(nargchk(2, 4, nargin));
error(nargoutchk(1, 1, nargout));

if nargin < 4 
    epsilon = 0;
end

if nargin < 3
    method = 'leave';
end

nd = ndims(img);
if  nd ~= 2  
    error('only 2D input supported.');
end

if ~isstruct(st)
    error('This is not a CTF structure.');
end

switch method
    case 'remove'
        method = 0;
    case 'leave'
        method = 1;
    case 'divide'
        method = 2;
    otherwise
        error('Method not supported.');
end

%cast the input variables to the correct type and call the xmipp function
img = double(img);
epsilon = double(epsilon);
method = int32(method);
img_out = tom_xmipp_ctf_correct_phase_wrapper(img,method,epsilon,st.K,st.Tm,st.kV,st.DeltafU,st.DeltafV,st.AzimuthalAngle,st.Cs,st.Ca,st.espr,st.ispr,st.alpha,st.DeltaF,st.DeltaR,st.Q0);