#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using hierarchical clustering principles
#
# Example use:
# ./xmipp_protocol_cl2d.py
#
# Author: Carlos Oscar Sanchez Sorzano, July 2011
#

# {begin_of_header}

# {eval} expandCommentRun()

#------------------------------------------------------------------------------------------------
# {section} Alignment parameters
#------------------------------------------------------------------------------------------------
# {file}{validate}(PathExists) Selfile with the input images:
""" This selfile points to the stack or metadata containing your images 
"""
InSelFile=''

# {file} Reference image:
""" The reference is optional. If no reference is provided the algorithm computes one of its own. """
ReferenceImage=""

# Maximum shift:
""" In pixels"""
MaxShift=10

# {expert} Number of iterations
""" Maximum number of iterations within each level
"""
NumberOfIterations=10

# {eval} expandParallel(threads=0)

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#		
# Main
#     
from protocol_cl2d_align import *

if __name__ == '__main__':
    protocolMain(ProtCL2DAlignment)
