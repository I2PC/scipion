#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Author: Carlos Oscar Oct 2013 
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
""" You may provide a very rough initial volume as a way to constraint the angular search.
    For instance, when reconstructing a fiber, you may provide a cylinder so that side views 
    are assigned to the correct tilt angle, although the rotational angle may be completely wrong
"""
InitialVolume = ''

#{expert} Number of simulated annealing iterations:
""" During the simulated annealing iterations, all those particles positively contributing to the improvement of the volume are considered.
In this way, the same image may participate several times from different projection directions (but different weights) depending
on whether it improves the correlation with the volume or not"""
NIterRandom = 10

#{expert} Initial temperature:
"""The initial temperature determines whether wrong orientations are considered or not. At the beginning of the simulated annealing
iterations it is important to take wrong directions in order to avoid local minima. You may increase this value and the number
of simulated annealing iterations as a way to have a slower cooling. """ 
T0 = 0.1

# {expert} Don't apply non-negative constraints
"""In between simulated annealing iterations, the reconstructed volume is constrained to have non-negative values. This helps
to reduce the search space, but it might be incorrect depending on the normalization of the input images"""
DontApplyPositiveConstraint=False

#{expert} Number of greedy iterations:
""" After simulating annealing, you may run some greedy iterations by simply applying a projection matching approach """
NIterGreedy = 2

#{expert} Percentage of rejected particles:
"""At each iteration, the lowest correlated particles are removed from the 3D reconstruction, although they may participate in the
next iteration""" 
Rejection = 50

#{expert} Angular sampling:
"""Angular sampling in degrees for generating the projection gallery""" 
AngularSampling = 5.0

# {eval} expandParallel(threads=0,mpi=0,hours=12)

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------

from protocol_initvolume_simanneal import *
#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtInitVolSimAnneal)
