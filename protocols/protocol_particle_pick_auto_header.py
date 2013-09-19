#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based automatic particle picking
#
# Author: Carlos Oscar Sorzano, September, 2011
#
#------------------------------------------------------------------------------------------------
# {begin_of_header} 

# {eval} expandCommentRun()

# {run}(particle_pick) Supervised particle picking run
"""
Select previous RUN of the <Supervised> particle picking.
"""
SupervisedRun = ""

# {file}(micrographs*.xmd) Set of micrographs to pick
""" Leave it empty if it is the same as in the supervised particle picking"""
MicrographsMd = ""

# {eval} expandJavaMemory()

# {eval} expandParallel(threads=0, hours=36)

#------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
#------------------------------------------------------------------------------------------------

from protocol_particle_pick_auto import *
#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtParticlePickingAuto)
