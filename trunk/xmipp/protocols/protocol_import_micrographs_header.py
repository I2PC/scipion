#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# {begin_of_header}

#{eval} expandCommentRun()

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

#------------------------------------------------------------------------------------------------
# {section}{has_question} Preprocess
#------------------------------------------------------------------------------------------------
# Do Preprocess
DoPreprocess=False

# Crop borders
""" Crop a given amount of pixels from each border.
    Set this option to -1 for not applying it."""
Crop = -1

# Remove bad pixels
""" Values will be thresholded to this multiple of standard deviations. Typical
    values are about 5, i.e., pixel values beyond 5 times the standard deviation will be
    substituted by the local median. Set this option to -1 for not applying it."""
Stddev = -1

# {wizard}(wizardBrowseCTF) Downsampling factor 
""" Set to 1 for no downsampling. Non-integer downsample factors are possible."""
Down = 1

#------------------------------------------------------------------------------------------------
# {section} Microscope description
#------------------------------------------------------------------------------------------------
# Microscope voltage (in kV)
Voltage = 200

# Spherical aberration (in mm)
SphericalAberration = 2.26

# Magnification rate
Magnification = 60000

# {list}(From image, From scanner) Sampling rate mode
SamplingRateMode="From image"

# {condition}(SamplingRateMode=="From image") Sampling rate (A/pixel)
SamplingRate = 1

# {condition}(SamplingRateMode=="From scanner") Scanned pixel size (in um/pixel)
ScannedPixelSize = 15

# {eval} expandParallel(threads=0,hours=6)

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_import_micrographs import *
if __name__ == '__main__':
    protocolMain(ProtImportMicrographs)
