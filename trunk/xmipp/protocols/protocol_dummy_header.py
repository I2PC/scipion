#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles, according to:
#
# Example use:
# ./xmipp_protocol_ml2d.py
#
#  Author:  Sjors Scheres, January 2008
# Updated:  J. M. de la Rosa Trevin July 2011
#
# {begin_of_header}
#{please_cite}
"""
for ML2D:  Scheres et al. (2005) J.Mol.Biol 348, 139-149
for MLF2D: Scheres et al. (2007) Structure 15, 1167-1177
"""
# {eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input
#-----------------------------------------------------------------------------
# {file}{view}{validate}(PathExists)Metadata to be splitted
""" Metadata to be splitted
"""
InputMd = "all_images.xmd"

# Randomize metadata before split?
""" Perform randomization of items before split
"""
Randomize = True

# {condition}(Randomize)
# Random seed
""" Just testing conditions
"""
RandomSeed = 3

# {validate}(IsInt) Number of part
""" Number of new metadatas splitted from the input one
"""
NumberOfParts = 2

# {validate}(NonEmpty) Part prefix
Prefix = "part"

# {eval} expandParallel(mpi=3, threads=1, condition="NumberOfParts>3")

#------------------------------------------------------------------------------------------------
# {section}{visualize} Visualization
#------------------------------------------------------------------------------------------------
# Visualize the class averages of all iterations in matrix-view?
DoMatrixAllIter=True
# Separately visualize class averages of the last iteration?
DoShowLastIter=True
# Plot model (and mirror) fractions of the last iteration?
DoShowFractions=True
# Plot convergence statistics for all iterations?
DoShowStatsAllIter=True


#------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
#------------------------------------------------------------------------------------------------

from protocol_dummy import *

if __name__ == '__main__':
    protocolMain(ProtDummy)
