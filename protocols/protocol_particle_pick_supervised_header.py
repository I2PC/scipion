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

# {run}(particle_pick) Manual particle picking RUN
""" 
Select previous RUN of the <Manual> particle picking.
"""
PickingRun = ""

# {wizard}(wizardChooseFamilyToExtractSupervised) Family manually picked
""" 
Select family within run
"""
Family = ""

# {eval} expandJavaMemory()

# Number of threads
""" 
This option provides shared-memory parallelization on multi-core machines.
It does not require any additional software, other than xmipp.
"""
NumberOfThreads = 2

# Fast picking
""" 
The fast version includes a Fourier filter while the non-fast version 
uses the Fourier filter and a Wavelet denoising. The fast version takes 
half the time of the non-fast version.
"""
Fast = True

# {expert} In-core picking
"""
If you can afford to load the micrographs in memory while picking, 
this option makes the program faster.
"""
InCore = False

#------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
#------------------------------------------------------------------------------------------------

from protocol_particle_pick_supervised import *
#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtParticlePickingSupervised)
