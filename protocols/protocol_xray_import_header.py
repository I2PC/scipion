#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of x-ray tomograms: 
#
# For each tomogram given, this script will perform 
# the requested operations below.
# For each tomogram a subdirectory will be created
#
# Author: Joaquin Oton, October 2013
#
#------------------------------------------------------------------------------------------------
# {begin_of_header}

#{eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input tomograms
#-----------------------------------------------------------------------------

# {eval}expandXrayImport()

# {eval} expandParallel(threads=1,mpi=0)

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_xray_import import *
if __name__ == '__main__':
    protocolMain(ProtXrayImport)
