#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# {begin_of_header}

#{eval} expandCommentRun()


#-----------------------------------------------------------------------------
# {section} Input micrographs
#-----------------------------------------------------------------------------
# {dir} Movies directory
"""Directory name from where to process all scanned movies"""
DirMovies = ''

# {wizard}(wizardMoviesExtension) Files to process
""" 
This is the movie files pattern, typically <*.mrcs> or <*.spi>
<Note:> you can use any wildcard like in shell expansion, e.g. <*3[1,2].mrcs>
"""
ExtMovies = ''

# {condition} Copy movies?
""" 
Set to <True> if you want to copy movies from input
folder to the protocol working directory
If you do preprocessing this option will be ignored.
"""
CopyMovies = False

# {condition} Average along each movie?
""" 
Set to <True> if you want to obtain an average 
for each movie.
"""
AverageMovies = True

#------------------------------------------------------------------------------------------------
# {section} Acquisition information
#------------------------------------------------------------------------------------------------
# Microscope voltage (in kV)
Voltage = 200

# Spherical aberration (in mm)
SphericalAberration = 2.26

# {list}(From image, From scanner) Sampling rate mode
SamplingRateMode = "From image"

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
from protocol_import_movies import *
if __name__ == '__main__':
    protocolMain(ProtImportMovies)
