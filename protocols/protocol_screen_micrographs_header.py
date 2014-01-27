#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of micrographs: 
#  - downsampling
#  - power spectral density (PSD) and CTF estimation on the micrograph
#
# For each micrograph given, this script will perform 
# the requested operations below.
# For each micrograph a subdirectory will be created
#
# Author: Carlos Oscar Sorzano, July 2011
#

# {begin_of_header}

#{eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} CTF Estimation
#-----------------------------------------------------------------------------
# {run}(import_micrographs,downsample_micrographs,emx_import_micrographs) Import Micrographs Run
ImportRun = ''

# {wizard}(wizardBrowseCTF2) CTF Downsampling factor 
""" 
Set to 1 for no downsampling. Non-integer downsample factors are possible.
This downsampling is only used for estimating the CTF and it does not affect 
any further calculation. Ideally the estimation of the CTF is optimal when 
the Thon rings are not too concentrated at the origin (too small to be seen)
and not occupying the whole power spectrum (since this downsampling might 
entail aliasing).
"""
DownsampleFactor = 1

# Automatic downsampling detection
""" If this option is chosen, the algorithm automatically tries by default the 
suggested Downsample factor; and if it fails, +1; and if it fails, -1.
"""
AutomaticDownsampling = True

# Amplitude Contrast
""" It should be a positive number, typically between 0.05 and 0.3."""
AmplitudeContrast = 0.1

# {wizard}(wizardBrowseCTF2) Lowest resolution
""" 
Give a value in digital frequency (i.e. between 0.0 and 0.5)
This cut-off prevents the typically peak at the center of the PSD 
to interfere with CTF estimation. The default value is 0.05, but for 
micrographs with a very fine sampling this may be lowered towards 0.0
"""
LowResolCutoff = 0.05

# {wizard}(wizardBrowseCTF2) Highest resolution
""" 
Give a value in digital frequency (i.e. between 0.0 and 0.5)
This cut-off prevents high-resolution terms where only noise exists 
to interfere with CTF estimation. The default value is 0.35, but it should
 be increased for micrographs with signals extending beyond this value.
However, if your micrographs extend further than 0.35, you should consider
sampling them at a finer rate.
"""
HighResolCutoff = 0.35

# {expert} Fast defocus estimate
FastDefocus = True

# {expert} Minimum defocus to search (in microns)
""" Minimum defocus value (in microns) to include in defocus search. 
Underfocus is represented by a positive number.
"""
MinFocus = 0.5

# {expert} Maximum defocus to search (in microns)
""" Maximum defocus value (in microns) to include in defocus search. 
Underfocus is represented by a positive number.
"""
MaxFocus = 10

# {expert} Window size
""" 
The PSD is estimated from small patches of this size. Bigger patches 
allow identifying more details. However, since there are fewer windows, 
estimations are noisier.
"""
WinSize = 256

# {expert} Automatic rejection
""" 
Reject micrographs meeting the following condition, e.g., ctfCritCorr13<0.4 OR ctfCritNormality<11. You may use
any of the following variables ctfDefocusU, ctfDefocusV, ctfCritFirstZero, ctfCritMaxFreq, ctfCritDamping, ctfCritfirstZeroRatio
ctfCritFitting, ctfCritCorr13, ctfCritPsdCorr90, ctfCritPsdInt, ctfCritPsdStdQ, ctfCritPsdPCA1, ctfCritPsdPCARuns, ctfCritNormality,
ctfCritFirstMinFirstZeroRatio, ctfCritCtfMargin, ctfCritNonAstigmaticValidty
"""
AutomaticRejection = "ctfCritFirstZero<5 OR ctfCritMaxFreq>20 OR ctfCritfirstZeroRatio<0.9 OR ctfCritfirstZeroRatio>1.1 OR ctfCritFirstMinFirstZeroRatio>10 OR ctfCritCorr13<0 OR ctfCritCtfMargin<2.5 OR ctfCritNonAstigmaticValidty<0.3"

# {eval} expandParallel(threads=0,hours=12)

# {eval} expandExpert()

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_screen_micrographs import *
if __name__ == '__main__':
    protocolMain(ProtScreenMicrographs)
