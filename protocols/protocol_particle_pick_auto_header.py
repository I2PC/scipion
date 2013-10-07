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

# {list}(Same as supervised, Other) Micrographs to pick
""" Select from which set of micrographs to pick using the training from supervised run.

If you use <Same as supervised>, the same set of micrographs used for training the picker
will be used at this point.

If you select <Other>, you can select another set of micrograph(normally from the same 
specimen) and pick them completely automatic using the trainned picker.
"""
MicrographSource = "Same as supervised"

# {condition}(MicrographSource=="Other"){file}(micrographs*.xmd) Set of micrographs to pick
""" Select other set of micrographs to pick using the trained picker"""
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
