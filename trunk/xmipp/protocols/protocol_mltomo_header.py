#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------

#
# Author: Carlos Oscar, December 2011
#
# {begin_of_header}

# {eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} General parameters
#-----------------------------------------------------------------------------
# {file} Metadata with volumes to align
VolumeList=''

# Generate references?
""" 
If you set to <No>, you should provide reference volumes (in a metadata file). 
"""
DoGenerateReferences = True

# {condition}(DoGenerateReferences) Number of references:
""" Number of references to be generated. """
NumberOfReferences = 3

# {file}{validate}(PathExists){condition}(not DoGenerateReferences) References metadata:
""" Metadata with the input reference volumes """
RefMd = "result_classes.xmd"

# {file} Missing regions:
""" Metadata with the description of the missing regions of the input volumes """
MissingMd = ""

#-----------------------------------------------------------------------------
# {section} Search parameters
#-----------------------------------------------------------------------------
# Angular sampling
""" Angular sampling in degrees for the different iterations """
AngSampling="15 7 4"

# Maximum angular change
""" Angular change from previous assignment """
AngLimit="360 21 8"

# Maximum shift change
""" Shift change from previous assignment """
ShiftLimit="-1 10 5"

#-----------------------------------------------------------------------------
# {section}{has_question} Advanced parameters
#-----------------------------------------------------------------------------
# Advanced parameters
SeeAdvanced = False

# Number of ML iterations to perform:
NumberOfIterations=25

# Symmetry:
""" See [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry]
    for a description of the symmetry groups format
    If no symmetry is present, give c1 """
Symmetry = "c1"

# Downsampled dimension
""" Internally, the volumes will be internally resized to this value (in pixels). Use a negative value to use the original volume size. 
"""
Dimension=-1

# Maximum resolution (in pixel^-1) to be used
""" The maximum (Nyquist) resolution is 0.5. Use smaller values, e.g. 0.45, to prevent high-resolution artifacts.
"""
MaximumResolution=0.35

# Apply small random perturbations on the angular sampling?
""" This option is recommended, as it will prevent the algorithm from getting stuck in local minima.
"""
DoPerturb=True

# Initial regularization parameter
""" We have found that imposing similarity on the difference references during the initial iterations of ML refinement
    may prevent the classification from getting stuck into local minima.
"""
InitialRegularization=0

# Number of iterations to use regularization
""" Useful values are 3-5. Note that during these iterations, the references will be enforced to be similar, so do not
   expect any useful classification until after these iterations. 
"""
NumberRegularizationSteps=5

# Additional xmipp_ml_tomo parameters:
""" For a complete description see the manual pages:
    http://http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Ml_tomo_v3
"""
ExtraParamsMLtomo=''

# {eval} expandParallel(threads=8, hours=24, mpi=2)
#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_mltomo import *
#        
# Main
#     
 
if __name__ == '__main__':
    protocolMain(ProtMLTomo)