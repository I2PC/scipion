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
""" Set of normal modes to explore 
"""
Modesfile=''

# Sampling rate (A/pix)
SamplingRate=2

#------------------------------------------------------------------------------------------------
# {section}{expert} Angular assignment and mode detection
#----------------------------------------------------------------------------------
# Trust region scale
"""This parameter scales the initial value of the trust region radius for optimization purposes.
Use larger values for larger expected deformation amplitudes."""
TrustRegionScale=1

# Use Projection Matching:
"""Use Projection Matching instead of Wavelets and Splines"""
ProjMatch=True

# Minimum angular sampling rate
MinAngularSampling=10

#{expert}{list_combo}(PCA, LTSA, DM, LLTSA, LPP, kPCA, pPCA, LE, HLLE) Projection method:
"""Projection methods: Principal Component Analysis, Local Tangent Space Alignment,
Diffusion map, Linear Local Tangent Space Alignment, Linearity Preserving Projection,
Kernel PCA, Probabilistic PCA, Laplacian Eigenmap, Hessian Locally Linear Embedding"""
ProjectionMethod="PCA"

#{expert} Output dimension:
"""Dimension of the lower dimensional space onto which the deformation coefficients are projected"""
OutputDim=2

# {eval} expandParallel(threads=0,mpi=2)

#------------------------------------------------------------------------------------------------
# {section}{visualize} Visualization
#------------------------------------------------------------------------------------------------

# Display raw deformation
""" Type 1 to see the histogram of raw deformation number 1; type 2 to see the histogram of raw deformation number 2, etc.
Type 1 2 to see the 2D plot of raw deformations number 1 vs 2. Type 1 2 3 to see the 3D plot of raw deformations 1, 2 and 3; etc.
""" 
DisplayRawDeformation="1"

# Display combined deformation
""" Type 1 to see the histogram of combined deformation number 1; type 2 to see the histogram of combined deformation number 2, etc.
Type 1 2 to see the 2D plot of combined deformations number 1 vs 2. Type 1 2 3 to see the 3D plot of combined deformations 1, 2 and 3; etc.
""" 
DisplayCombinedDeformation="1"

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#		
# Main
#     
from protocol_nma_alignment import *

if __name__ == '__main__':
    protocolMain(ProtNMAAlignment)
