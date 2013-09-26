#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of single-particles: 
# Author: Carlos Oscar, August 2011
#
# {begin_of_header}

# {eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input
#-----------------------------------------------------------------------------
# {file}(images.xmd){validate}(PathExists) Input images:
""" 
This selfile points to the stack or metadata containing your images 
"""
InSelFile = ''

# {eval} expandThreshold()

# {eval} expandFilter()

# {eval} expandMask()

# {eval} expandResize()

# {eval} expandParallel(mpi=8, threads=0, hours=6)
#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_preprocess_particles import *
#        
# Main
#     
 
if __name__ == '__main__':
    # create preprocess_particles_class object
    protocolMain(ProtPreprocessParticles)
