#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for flexible angular alignment
#
# Author: Carlos Oscar Sanchez Sorzano, May 2013
#         Qiyu Jin
#         Slavica Jonic
#

# {begin_of_header}

# {eval} expandCommentRun()

#------------------------------------------------------------------------------------------------
# {section} General parameters
#------------------------------------------------------------------------------------------------
# {file}(images*.xmd){validate}(PathExists) Selfile with the input images:
""" This selfile points to the stack or metadata containing your images 
"""
InSelFile=''

# {file}(*.pdb){validate}(PathExists) PDB to apply the modes to:
""" Atomic or pseudo-atomic structure to apply the normal modes 
"""
PDBfile=''

# {file}(modes*.xmd){validate}(PathExists) Normal modes:
""" Set of normal modes to explore. The Normal Mode Analysis protocol can be used to obtain this file in the required format. 
"""
Modesfile=''

# Sampling rate (A/pix)
SamplingRate=2

#------------------------------------------------------------------------------------------------
# {section}{expert} Angular assignment and mode detection
#----------------------------------------------------------------------------------
# Trust region scale
"""This parameter can be used to scale the default initial value of the trust-region radius in the elastic-alignment optimizer, if much larger movements are expected than those that are usually encountered. The default value of 1 works well for usually encountered movement amplitudes."""
TrustRegionScale=1

# Use Projection Matching:
"""If the answer to this question is "No" (default), Wavelets-and-Splines strategy will be used for rigid-body alignment. If the answer is "Yes", Projection Matching will be used for rigid-body alignment. Preferably use Wavelets-and-Splines (more accurate) instead of Projection Matching (faster). In the case of Wavelets-and-Splines, the size of images should be a power of 2."""
ProjMatch=False

# Discrete angular sampling rate
"""This parameter is used in Projection Matching and Wavelets methods for a rough rigid-body alignment. It is the angular step (in degrees) with which the library of reference projections is computed. This alignment is refined with Splines method if Wavelets-and-Splines alignment is chosen."""
DiscreteAngularSampling=10

# {eval} expandParallel(threads=0,mpi=2)

#------------------------------------------------------------------------------------------------
# {section}{visualize} Visualization
#------------------------------------------------------------------------------------------------

# Display raw deformation
""" Type 1 to see the histogram of raw deformation number 1; type 2 to see the histogram of raw deformation number 2, etc.
Type 1 2 to see the 2D plot of raw deformations number 1 vs 2. Type 1 2 3 to see the 3D plot of raw deformations 1, 2 and 3; etc.
""" 
DisplayRawDeformation="1"

#{view} Analyze with MATLAB
""" In MATLAB you may create image clusters and draw deformation trajectories
""" 
AnalyzeMATLAB=True

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#		
# Main
#     
from protocol_nma_alignment import *

if __name__ == '__main__':
    protocolMain(ProtNMAAlignment)
