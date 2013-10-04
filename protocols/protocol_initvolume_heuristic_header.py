#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Author: Carlos Oscar Sorzano, March 2013 
#
# {begin_of_header}

#{eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input
#-----------------------------------------------------------------------------

# {file}(classes*.xmd){validate}(PathExists) Set of classes:
""" 
Provide a metadata or stack with classes
"""
Classes = ""

# Symmetry group
""" See [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry]
    for a description of the symmetry groups format
    If no symmetry is present, give c1
"""
SymmetryGroup ='c1'

#{file}(*.vol) Initial Volume:
""" Initial volume to be used 
"""
InitialVolume = ''

#{expert} Number of heuristic iterations
""" During the heuristic iterations, all those particles positively contributing to the improvement of the volume are considered.
In this way, the same image may participate several times from different projection directions (but different weights) depending
on whether it improves the correlation with the volume or not"""
NIterRandom = 7

#{expert} Number of greedy iterations
""" These are standard projection matching iterations performed after the heuristic iterations"""
NIterGreedy = 3

#{expert} Percentage of rejected particles
"""At each iteration, the lowest correlated particles are removed from the 3D reconstruction, although they may participate in the
next iteration""" 
Rejection = 25

#{expert} Apply positive constraint
"""If the classes are properly normalized (it suffices that you normalized the extracted images before calculating classes), then
adding the constraint that the volume should have positive values normally helps to have better reconstructions.""" 
Positive = True

#{expert} Keep intermediate volumes
"""Keep all the volumes of the different iterations""" 
KeepIntermediate = False

# {eval} expandParallel(threads=4,mpi=0,hours=12)

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------

from protocol_initvolume_heuristic import *
#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtInitVolH)
