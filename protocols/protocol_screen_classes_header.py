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

# {file}(volume*.vol){validate}(PathExists) Volume to compare classes to:
""" 
Provide a volume against which classes will be compared
"""
Volume = ""

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

# {eval} expandParallel(threads=0,hours=12)

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------

from protocol_screen_classes import *
#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtScreenClasses)
