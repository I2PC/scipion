#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for creaing a symmetric initial volume
#
# Example use:
# ./xmipp_protocol_rct.py
#
# Author: Carlos Oscar Sorzano, March 2013 
#
# {begin_of_header}

#{eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input
#-----------------------------------------------------------------------------

# {file}{validate}(PathExists) Side views:
""" 
Provide a metadata or stack with side views
"""
SideViews = ""

# {file} Top views:
""" 
Provide a metadata or stack with top views. This field is optional
"""
TopViews = ""

# Symmetry group
""" See [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry]
    for a description of the symmetry groups format
    If no symmetry is present, give c1
"""
SymmetryGroup ='c6'

#{expert} Number of volumes
""" Number of volumes to construct """
Nvolumes=20

# {eval} expandParallel(threads=0,hours=12)

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------

from protocol_symmetric_initial import *
#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtSymmetric)
