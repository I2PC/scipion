#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for building initial references by common lines
# This protocol is based on EMAN 1
#
# Example use:
# ./xmipp_protocol_commonlines.py
#
# Author:Carlos Oscar Sorzano, January 2011
#

# {begin_of_header}

# {eval} expandCommentRun()

#------------------------------------------------------------------------------------------------
# {section} Common line parameters
#------------------------------------------------------------------------------------------------
# {file}{validate}(PathExists) Selfile with the input class averages:
""" This selfile points to the stack or metadata containing your images 
"""
InSelFile=''

# Particle radius
Radius=25
""" Maximum radius of the particle
"""
# Symmetry
""" c1=No symmetry; c2, c3, ...
"""
Symmetry='c1'

#{eval} expandParallel(threads=0)

#-----------------------------------------------------------------------------
# {end_of_header} do not change anything bellow this line unless you know what you are doing
#-----------------------------------------------------------------------------

from protocol_commonlines import *
  
if __name__ == '__main__':
    protocolMain(ProtCommonLines)
