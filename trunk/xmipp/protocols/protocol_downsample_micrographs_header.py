#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# {begin_of_header}

#{eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input micrographs
#-----------------------------------------------------------------------------
# {run}(import_micrographs) Import Micrographs Run
ImportRun = ''

# {wizard}(wizardBrowseCTF) Downsampling factor 
""" Non-integer downsample factors are possible. Must be larger than 1. """
DownsampleFactor = 2

# {eval} expandParallel(threads=0,hours=6)

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_downsample_micrographs import *
if __name__ == '__main__':
    protocolMain(ProtDownsampleMicrographs)
