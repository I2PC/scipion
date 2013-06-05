#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for screening a number of classes comparing them to a volume
#
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

#{expert} Angular sampling
""" In degrees. This sampling defines how fine the projection gallery from the volume is explored.
"""
AngularSampling = 5

#{expert} Number of RANSAC iterations
""" Number of initial volumes to test by RANSAC
"""
NRansac = 300

#{expert} Number of grids per dimension 
""" Number of grids per dimension. The total number of classes used consists in n x n
"""
NumGrids = 3

#{expert} Inlier percentil
""" correlation value using the estimated initial volume to determine if an experimental projection is an
inlier or outlier
"""
CorrThresh = 0.8

#{expert} Number of best volumes to refine
""" Number of best volumes to refine 
"""
NumVolumes =3

#{expert} Number of iterations to perform to refine the volumes
""" Number of iterations to perform to refine the volumes 
"""
NumIter = 10

#{file}(*.vol) Initial Volume:
""" Initial volume to be used 
"""
InitialVolume = ''

#{expert} Max frequency of the initial volume
""" Max frequency of the initial volume in Angstroms
"""
MaxFreq = 40

# Sampling Rate
""" Sampling rate (A/px)
"""
Ts = '1'

#{expert} Use all images to refine
""" When refining a RANSAC volume, use all images to refine it instead of only inliers
"""
UseAll=False

# {eval} expandParallel(threads=0,hours=12)

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------

from protocol_initvolume_ransac import *
#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtInitVolRANSAC)
