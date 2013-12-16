#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------

#
# Author: Carlos Oscar, September 2013
#
# {begin_of_header}

# {eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} General parameters
#-----------------------------------------------------------------------------
# {file}(*.xmd) Metadata with volumes to align
VolumeList=''

# Number of classes:
""" Number of references to be generated. """
NumberOfReferences = 3

# {expert} Number of CL iterations to perform:
NumberOfIterations=15

# Generate aligned volumes
GenerateAligned=True

# Don't align
"""Volumes are already aligned, only classify"""
DontAlign=False

#-----------------------------------------------------------------------------
# {section}{has_question} Initial classes
#-----------------------------------------------------------------------------
# See parameters
SeeInitial = False

# Let CLTomo generate initial classes
DoGenerateInitial = True

# {condition}(DoGenerateInitial) Initial number of classes:
""" The algorithm proceeds by doubling the initial number of classes till the final number of classes is reached."""
NumberOfReferences0 = 1

# {condition}(DoGenerateInitial) Randomize initial volume orientation:
""" Use this option if all the input volumes have the same missing wedge or if they have not been previously aligned"""
RandomizeOrientation = False

# {file}(classes*.xmd){validate}(PathExists){condition}(not DoGenerateInitial) Initial classes:
""" Metadata with the input reference volumes """
RefMd = ""

#-----------------------------------------------------------------------------
# {section}{has_question} Search limits
#-----------------------------------------------------------------------------
# See limits
SeeLimits = False

# Maximum rotational angle
MaxRot=360

# Maximum tilt angle
MaxTilt=360

# Maximum in-plane angle
MaxPsi=360

# Maximum shift in X
MaxShiftX=5

# Maximum shift in Y
MaxShiftY=5

# Maximum shift in Z
MaxShiftZ=5

#-----------------------------------------------------------------------------
# {section}{has_question} Constraints
#-----------------------------------------------------------------------------
# See constraints
SeeConstraints = False

# Symmetry:
""" See [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry]
    for a description of the symmetry groups format
    If no symmetry is present, give c1 """
Symmetry = "c1"

# {file}(mask*.vol){validate}(PathExists) Mask:
""" Provide any mask parameters valid after --mask see [http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Transform_mask_v3]"""
Mask = ""

# Maximum resolution (in pixel^-1) to be used
""" The maximum (Nyquist) resolution is 0.5. Use smaller values, e.g. 0.45, to prevent high-resolution artifacts.
"""
MaximumResolution=0.25

# Percentage of Fourier coefficients to drop
""" A value of 97.5 drops 97.5% of the smallest Fourier coefficients
"""
Sparsity=97.5

# Percentage of Wavelet coefficients to drop
""" A value of 97.5 drops 97.5% of the smallest Wavelet coefficients
"""
DWTSparsity=99.0

# {eval} expandParallel(threads=0, hours=24, mpi=8)
#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_cltomo import *
#        
# Main
#     
 
if __name__ == '__main__':
    protocolMain(ProtCLTomo)
