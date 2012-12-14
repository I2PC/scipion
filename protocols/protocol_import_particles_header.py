#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based importing single-particles: 
#
# Author: Carlos Oscar, August 2011
#
# {begin_of_header}

# {eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input
#-----------------------------------------------------------------------------
# {file}{validate}(PathExists) Input images:
InputFile = ''

# {validate}(IsFloat) Sampling rate (A/pixel)
SamplingRate = 1

# Copy images?
""" 
If set to Yes, the images will be copied from their original location 
This option is ignored if you do preprocess"""
DoCopy = False

# Import all images from input?
ImportAll = True

# {list_combo}(First particles,Random particles){condition}(not ImportAll) Subset:
SubsetMode="First particles"

# {condition}(not ImportAll) Number of particles:
""" A subset with the first N particles will be imported """
Nsubset = 500

# Sort particles?
"""
If set to Yes, the particles will be sorted by similarity.
A ZScore column will be added to output <images.xmd>.
A higher ZScore means a more different particles.
""" 
DoSort = False

# {eval} expandPreprocessFilterMask(allowFlip=False)

# {eval} expandParallel(threads=0, mpi=8, hours=12, condition="ImportAll")

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_import_particles import *
#        
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtImportParticles)