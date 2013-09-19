#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# {begin_of_header}

#{eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input micrographs
#-----------------------------------------------------------------------------
# {run}(import_micrographs,screen_micrographs,emx_import) Micrographs to downsample
""" List with input micrographs. It is obtained from an execution of
either <import_micrograph> or <screen_micrograph>. <BE AWARE>, if you
use as input an execution of <import_micrographs> there will be not
<CTF> information available."""
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
