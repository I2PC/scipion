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
#------------------------------------------------------------------------------------------
# {section}{has_question} Comment
#------------------------------------------------------------------------------------------
# Display comment
DisplayComment = False

# {text} Write a comment:
""" 
Describe your run here...
"""
#-----------------------------------------------------------------------------
# {section} Run parameters
#-----------------------------------------------------------------------------
# Run name:
""" This will identify your protocol run. It need to be unique for each protocol. You could have run1, run2 for protocol X, but not two
run1 for it. This name together with the protocol output folder will determine the working dir for this run.
"""
RunName = "run_001"

# {list}(Resume, Restart) Run behavior
""" Resume from the last step, restart the whole process
"""
Behavior = "Resume"

#-----------------------------------------------------------------------------
# {section} Input micrographs
#-----------------------------------------------------------------------------
# {dir} Micrographs directory
"""Directory name from where to process all scanned micrographs"""
DirMicrographs = 'Micrographs'

# Files to process
""" This is typically *.tif or *.ser, but may also be *.mrc, *.spi 
    Note that any wildcard is possible, e.g. *3[1,2].tif
"""
ExtMicrographs = '*.mrc'

# Do preprocess
DoPreprocess = False

# Perform CTF estimation
DoCtfEstimate=True

#------------------------------------------------------------------------------------------------
# {section}{condition}(DoPreprocess) Preprocess
#------------------------------------------------------------------------------------------------
# Crop borders
""" Crop a given amount of pixels from each border.
    Set this option to -1 for not applying it."""
Crop = -1

# Remove bad pixels
""" Values will be thresholded to this multiple of standard deviations. Typical
    values are about 5, i.e., pixel values beyond 5 times the standard deviation will be
    substituted by the local median. Set this option to -1 for not applying it."""
Stddev = -1

# {wizard}(wizardBrowseJCTF) Downsampling factor 
""" Set to 1 for no downsampling. Non-integer downsample factors are possible."""
Down = 1

#------------------------------------------------------------------------------------------------
# {section}{condition}(DoCtfEstimate) CTF estimation
#------------------------------------------------------------------------------------------------

# Microscope voltage (in kV)
Voltage = 200

# Spherical aberration (in mm)
SphericalAberration = 2.26

# Magnification rate
Magnification = 70754

# {list}(From image, From scanner) Sampling rate mode
SamplingRateMode="From image"

# {condition}(SamplingRateMode=="From image") Sampling rate (A/pixel)
SamplingRate = 1

# {condition}(SamplingRateMode=="From scanner") Scanned pixel size (in um/pixel)
ScannedPixelSize = 15

# Amplitude Contrast
""" It should be a negative number, typically between -0.05 and -0.3."""
AmplitudeContrast = -0.1

# {wizard}(wizardBrowseJCTFMeasure){expert} Lowest resolution for CTF estimation
""" Give a value in digital frequency (i.e. between 0.0 and 0.5)
    This cut-off prevents the typically peak at the center of the PSD to interfere with CTF estimation.  
    The default value is 0.05, but for micrographs with a very fine sampling this may be lowered towards 0.0
"""
LowResolCutoff = 0.05

# {wizard}(wizardBrowseJCTFMeasure){expert} Highest resolution for CTF estimation
""" Give a value in digital frequency (i.e. between 0.0 and 0.5)
    This cut-off prevents high-resolution terms where only noise exists to interfere with CTF estimation.  
    The default value is 0.35, but it should be increased for micrographs with signals extending beyond this value.
    However, if your micrographs extend further than 0.35, you should consider sampling them at a finer rate.
"""
HighResolCutoff = 0.35

# {expert} Minimum defocus to search (in microns)
""" Minimum defocus value (in microns) to include in defocus search
"""
MinFocus = 0.5

# {expert} Maximum defocus to search (in microns)
""" Maximum defocus value (in microns) to include in defocus search
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

#------------------------------------------------------------------------------------------
# {section}{condition}(DoCtfEstimate or DoPreprocess) Parallelization
#------------------------------------------------------------------------------------------
# Number of MPI processes
""" Set to 1 if you do not have MPI installed"""
NumberOfMpi = 8

# Submit to queue
"""Submit to queue
"""
SubmitToQueue = False

# {condition}(SubmitToQueue) Queue name
"""Name of the queue to submit the job
"""
QueueName = "default"

# {condition}(SubmitToQueue) Queue hours
"""This establish a maximum number of hours the job will
be running, after that time it will be killed by the
queue system
"""
QueueHours = 6

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
#        
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtScreenMicrographs)
