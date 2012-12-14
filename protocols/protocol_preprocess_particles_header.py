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
# {file}(images*.xmd){validate}(PathExists) Selfile with the input images:
""" This selfile points to the stack or metadata containing your images 
"""
InSelFile = ''

# Scale
""" Change the scale of the input images """
DoScale = False

#{condition}(DoScale) New image size
NewSize = 0

# Crop
"""
This is the desired output size(in pixels) after cropping.
"""
DoCrop = False

# {condition}(DoCrop) Output size:
""" 
In pixels
"""
CropSize = 64

# {eval} expandParticlesPreprocess(allowFlip=False)

# {eval} expandFilter()

# {eval} expandMask()

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
