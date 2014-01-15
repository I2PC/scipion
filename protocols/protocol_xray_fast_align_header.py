#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based fast alignment of x-ray tomograms using IMOD.
#
# IMOD steps provided by:
#                         Javier Chichon
#                         Javier Conesa
#
# Author: Joaquin Oton, October 2013
#
#------------------------------------------------------------------------------------------------
# {begin_of_header}

#{eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input tomograms
#-----------------------------------------------------------------------------
# {run}(xray_import) Import Tomograms Run
ImportRun = ''

# {eval}expandXrayFastAlign()

# {eval} expandParallel(threads=0, mpi=0)

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_xray_fast_align import *
if __name__ == '__main__':
    protocolMain(ProtXrayFastAlign)
