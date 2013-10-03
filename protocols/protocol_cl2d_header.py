#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using hierarchical clustering principles
#
# Example use:
# ./xmipp_protocol_cl2d.py
#
# Author: Carlos Oscar Sanchez Sorzano, July 2009
#

# {begin_of_header}

# {eval} expandCommentRun()

#------------------------------------------------------------------------------------------------
# {section} CL2D parameters
#------------------------------------------------------------------------------------------------
# {file}(images*.xmd){validate}(PathExists) Selfile with the input images:
""" This selfile points to the stack or metadata containing your images 
"""
InSelFile=''

# {wizard}(wizardCL2DNumberOfClasses) Number of references (or classes) to be used:
NumberOfReferences=64

# {expert} Number of initial references
""" Initial number of initial models
"""
NumberOfInitialReferences=4

# {expert} Number of iterations
""" Maximum number of iterations within each level
"""
NumberOfIterations=10

# {expert}{list_combo}(correlation, correntropy) Comparison method
""" Use correlation or correntropy """
ComparisonMethod='correlation'

# {expert}{list_combo}(classical, robust) Clustering criterion
""" Use the classical clustering criterion or the robust clustering criterion """
ClusteringMethod='classical'

# {expert} Additional parameters for classify_CL2D
""" --verbose, --corrSplit, --dontMirrorImages ... see [http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Classify_mpi_cl2d_v3]
"""
AdditionalParameters=''

#------------------------------------------------------------------------------------------------
# {section}{expert}{condition}(NumberOfReferences>NumberOfInitialReferences) Core analysis
#------------------------------------------------------------------------------------------------
# Junk Zscore
""" Which is the average Z-score to be considered as junk. Typical values
    go from 1.5 to 3. For the Gaussian distribution 99.5% of the data is
    within a Z-score of 3. Lower Z-scores reject more images. Higher Z-scores
    accept more images."""
thZscore=3

# PCA Zscore
""" Which is the PCA Z-score to be considered as junk. Typical values
    go from 1.5 to 3. For the Gaussian distribution 99.5% of the data is
    within a Z-score of 3. Lower Z-scores reject more images. Higher Z-scores
    accept more images.
    
    This Z-score is measured after projecting onto the PCA space."""
thPCAZscore=3

# {validate}(IsInt) Tolerance
""" An image belongs to the stable core if it has been with other images in the same class
in all the previous levels except possibly a few of them. Tolerance defines how few is few.
Tolerance=0 means that an image must be in all previous levels with the rest of images in
the core."""
Tolerance=1

# {eval} expandParallel(threads=0)

#------------------------------------------------------------------------------------------------
# {section}{visualize} Visualization
#------------------------------------------------------------------------------------------------
# {list_combo}(Classes,Class Cores,Class Stable Cores) What to show
WhatToShow="Class Cores"

# Visualize class hierarchy
DoShowHierarchy=True

# Visualize last level
DoShowLast=True

# {condition}(DoShowLast==False) Visualize several levels
"""Create a  list of levels like: 0,1,3 or 0-3 """
LevelsToShow=""

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#		
# Main
#     
from protocol_cl2d import *

if __name__ == '__main__':
    protocolMain(ProtCL2D)
