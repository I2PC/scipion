#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based automatic particle picking
#
# Author: Carlos Oscar Sorzano, September, 2011
#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {begin_of_header} 

# {eval} expandCommentRun()

#------------------------------------------------------------------------------------------
# {section} Picking parameters
#------------------------------------------------------------------------------------------
# {run}(particle_pick) Manual particle picking run
""" Directory with the manual picking
"""
PickingRun = ""

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
