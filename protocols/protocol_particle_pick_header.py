#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based manual particle picking
#
# For each micrograph in the PreprocessingDir, this program will launch
# the xmipp_mark program 
# A graphical interface exists to identify micrographs that have been finished
#
# Author: Sjors Scheres, March 2007
# Author: Carlos Oscar Sorzano, June, 2011
#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {begin_of_header} 

# {eval} expandCommentRun()

# {run}(downsample_micrographs,screen_micrographs,import_micrographs) Micrographs RUN
"""
Select desired RUN from which you obtained the micrographs.

Possible input protocols are:
<Import Micrographs>
<Screen Micrographs> or
<Downsample Micrographs>

If you want to correct for the CTF in a subsequent step, you must use a screen or downsample run.
"""
ImportRun = ""


# {hidden} launchGUI
LaunchGUI = True



#------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
#------------------------------------------------------------------------------------------------

from protocol_particle_pick import *
#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtParticlePicking)
