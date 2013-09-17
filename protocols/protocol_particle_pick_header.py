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

# {eval} expandJavaMemory()

# {expert} Number of threads
""" 
This option provides shared-memory parallelization on multi-core machines.
It does not require any additional software, other than xmipp.
"""
NumberOfThreads = 2

# {expert} Fast picking
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

from protocol_particle_pick import *
#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtParticlePicking)
