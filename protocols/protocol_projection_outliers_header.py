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

# {file}(images*.xmd){validate}(PathExists) Set of images:
""" 
Provide a metadata or stack with images
"""
Images = ""

# {file}(volume*.vol){validate}(PathExists) Volume to compare images to:
""" 
Provide a volume against which images will be compared
"""
Volume = ""

# Volume has been corrected by CTF:
VolumeIsCTFCorrected = False

# Symmetry group
""" See [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry]
    for a description of the symmetry groups format
    If no symmetry is present, give c1
"""
SymmetryGroup ='c1'

#{expert} Angular sampling
""" In degrees. This sampling defines how fine the projection gallery from the volume is explored.
"""
AngularSampling = 3

# {eval} expandParallel(threads=0,hours=12)

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------

from protocol_projection_outliers import *
#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtProjectionOutliers)
