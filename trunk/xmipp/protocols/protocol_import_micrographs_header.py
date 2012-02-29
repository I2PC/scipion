#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# {begin_of_header}

#{eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input micrographs
#-----------------------------------------------------------------------------
# {dir} Micrographs directory
"""Directory name from where to process all scanned micrographs"""
DirMicrographs = 'InputData'

# Files to process
""" 
This is the micrographs files pattern, typically <*.tif> or <*.ser>, 
but may also be <*.mrc> or <*.spi> 
<Note:> you can use any wildcard like in shell expansion, e.g. <*3[1,2].tif>
"""
ExtMicrographs = '*.mrc'

# Are micrographs tilt pairs?
TiltPairs = False

# {condition}(TiltPairs){wizard}(wizardTiltPairs) Pair assignment
""" 
You should provide a metadata with entries for each 
untilted and tilted pairs. It should reference to micrographs
in previous selected folder. The wizard will help you
to create this metadata.
"""
PairDescr = ""

# Preprocess micrographs?
"""Perform some preprocessing operations on micrographs"""
DoPreprocess = False

# {condition}(not DoPreprocess) Copy micrographs?
""" 
Set to <True> if you want to copy micrographs from input
folder to the protocol working directory
"""
CopyMicrographs = False

#------------------------------------------------------------------------------------------------
# {condition}(DoPreprocess){section} Preprocessing
#------------------------------------------------------------------------------------------------

# Crop borders?
""" Crop a given amount of pixels from each border. """
DoCrop = False

# {condition}(DoCrop) Pixels to crop
""" Amount of pixels you want to crop from borders """
Crop = 10

# Remove bad pixels?
""" 
Values will be thresholded to this multiple of standard deviations. 
Typical values are about 5, i.e., pixel values beyond 5 times the 
standard deviation will be substituted by the local median. 
Set this option to <-1> for not applying it."""
DoRemoveBadPixels = False

# {condition} (DoRemoveBadPixels) Multiple of Stddev
Stddev = 5

# Downsample micrographs?
DoDownsample = False
# {condition}(DoDownsample){wizard}(wizardBrowseCTF) Downsampling factor 
""" Set to 1 for no downsampling. Non-integer downsample factors are possible."""
DownsampleFactor = 1

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
