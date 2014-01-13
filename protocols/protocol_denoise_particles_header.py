#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using hierarchical clustering principles
#
# Example use:
# ./xmipp_protocol_cl2d.py
#
# Author: Carlos Oscar Sanchez Sorzano, July 2011
#

# {begin_of_header}

# {eval} expandCommentRun()

#------------------------------------------------------------------------------------------------
# {section} Basis construction
#------------------------------------------------------------------------------------------------
# {file}(classes*.xmd){validate}(PathExists) Input classes:
""" This selfile points to the stack or metadata containing your classes. Make sure you point
    to the block containing the classes
"""
InClasses=''

# {expert} Max. number of classes
""" Maximum number of classes
"""
MaxNumberOfClasses=128

# {expert} Number of PCA bases
""" Number of PCA bases
"""
NumberOfBases=200

#------------------------------------------------------------------------------------------------
# {section} Denoising
#------------------------------------------------------------------------------------------------
# {file}(images*.xmd){validate}(PathExists) Input images:
""" This selfile points to the stack or metadata containing your images. It is important that
    the images have alignment information with respect to the chosen set of classes.
    This is the standard situation after CL2D or ML2D.
"""
InStack=''

# {expert} Number of PCA bases on which to project
""" Number of PCA bases
"""
NumberOfBasesDenoising=200

# {eval} expandParallel(threads=0)

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#		
# Main
#     
from protocol_denoise_particles import *

if __name__ == '__main__':
    protocolMain(ProtDenoiseParticles)
