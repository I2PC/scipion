#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# {begin_of_header}

#{eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input micrographs
#-----------------------------------------------------------------------------
# {dir} Micrographs directory
"""Directory name from where to process all scanned micrographs"""
DirMicrographs = ''

# {wizard}(wizardMicrographExtension) Files to process
""" 
This is the micrographs files pattern, typically <*.tif> or <*.ser>, 
but may also be <*.mrc> or <*.spi> 
<Note:> you can use any wildcard like in shell expansion, e.g. <*3[1,2].tif>
"""
ExtMicrographs = ''

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

# {expert} Check image corners for problems
"""After doing all the preprocess check that image borders have similar variance"""
DoCheckBorders = True

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

# Take Logarithm ?
""" 
Depending on your acquisition system you may have to take the logarithm
or not in order to have a linear relationship between the gray values
in the image and those in the volume. (log (x+1)) """
DoLog = False 

# Remove bad pixels?
""" 
Values will be thresholded to this multiple of standard deviations. 
Typical values are about 5, i.e., pixel values beyond 5 times the 
standard deviation will be substituted by the local median. 
Set this option to <-1> for not applying it."""
DoRemoveBadPixels = False

# {condition} (DoRemoveBadPixels) Multiple of Stddev
Stddev = 5

#------------------------------------------------------------------------------------------------
# {section} Microscope description
#------------------------------------------------------------------------------------------------
# Microscope voltage (in kV)
Voltage = 200

# Spherical aberration (in mm)
SphericalAberration = 2.26

# {list}(From image, From scanner) Sampling rate mode
SamplingRateMode="From image"

# {condition}(SamplingRateMode=="From image"){validate}(IsFloat) Sampling rate (A/pixel)
SamplingRate = ""

# {condition}(SamplingRateMode=="From scanner") Magnification rate
Magnification = 60000
# {condition}(SamplingRateMode=="From scanner"){validate}(IsFloat) Scanned pixel size (in um/pixel)
ScannedPixelSize = ""

# {eval} expandParallel(threads=0,hours=6)

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_import_micrographs import *
if __name__ == '__main__':
    protocolMain(ProtImportMicrographs)
