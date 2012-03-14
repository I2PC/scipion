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
# {run}(import_micrographs) Import Micrographs Run
ImportRun = ''

# {wizard}(wizardBrowseCTF2) Downsampling factor 
""" Set to 1 for no downsampling. Non-integer downsample factors are possible."""
DownsampleFactor = 1

# Amplitude Contrast
""" It should be a positive number, typically between 0.05 and 0.3."""
AmplitudeContrast = 0.1

# {wizard}(wizardBrowseCTF2) Lowest resolution for CTF estimation
""" Give a value in digital frequency (i.e. between 0.0 and 0.5)
    This cut-off prevents the typically peak at the center of the PSD to interfere with CTF estimation.  
    The default value is 0.05, but for micrographs with a very fine sampling this may be lowered towards 0.0
"""
LowResolCutoff = 0.05

# {wizard}(wizardBrowseCTF2) Highest resolution for CTF estimation
""" Give a value in digital frequency (i.e. between 0.0 and 0.5)
    This cut-off prevents high-resolution terms where only noise exists to interfere with CTF estimation.  
    The default value is 0.35, but it should be increased for micrographs with signals extending beyond this value.
    However, if your micrographs extend further than 0.35, you should consider sampling them at a finer rate.
"""
HighResolCutoff = 0.35

# {expert} Minimum defocus to search (in microns)
""" Minimum defocus value (in microns) to include in defocus search. Underfocus is represented by a positive number.
"""
MinFocus = 0.5

# {expert} Maximum defocus to search (in microns)
""" Maximum defocus value (in microns) to include in defocus search. Underfocus is represented by a positive number.
"""
MaxFocus = 10

# {expert} Window size
""" The PSD is estimated from small patches of this size. Bigger patches allow identifying more details.
    However, since there are fewer windows, estimations are noisier.
"""
WinSize = 256

# Do CTFFIND
""" This option can be applied if ctffind3.exe is on the system PATH.
"""
DoCtffind = False

# {condition}(DoCtffind){expert} Defocus step for CTFFIND (in microns)
""" Step size for defocus search (in microns)
"""
StepFocus = 0.1

# {eval} expandParallel(threads=0,hours=12)

# {hidden} Show expert options
"""If True, expert options will be displayed
"""
ShowExpertOptions = False

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_screen_micrographs import *
if __name__ == '__main__':
    protocolMain(ProtScreenMicrographs)
