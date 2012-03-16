#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for image classification using self-organizing maps
#
# Example use:
# ./xmipp_protocol_kerdensom.py
#
# Author:Carlos Oscar Sorzano, January 2011
#
#
# {begin_of_header}

# {eval} expandCommentRun()

#------------------------------------------------------------------------------------------------
# {section} KerdenSOM parameters
#------------------------------------------------------------------------------------------------
# {file}{validate}(PathExists) Selfile with the input images:
""" This selfile points to the stack or metadata containing your images 
"""
InSelFile=''

# Use mask
UseMask=False

# {file}{condition}(UseMask) Mask file 
Mask=''

# X-dimension of the map:
SomXdim=7

# Y-dimension of the map:
SomYdim=7

# {expert} Initial regularization factor:
""" The kerdenSOM algorithm anneals from an initial high regularization factor
    to a final lower one, in a user-defined number of steps.
    If the output map is too smooth, lower the regularization factors
    If the output map is not organized, higher the regularization factors
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM
"""
SomReg0=1000

# {expert} Final regularization factor:
SomReg1=200

# {expert} Regularization steps
"""Number of steps to lower the regularization factor"""
SomSteps=5

# {expert} Additional kerdenSOM parameters:
""" For a complete description 
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM
"""
KerdensomExtraCommand=''

# {eval} expandParallel(threads=0, mpi=0, hours=96)

#------------------------------------------------------------------------------------------------
# {end_of_header} do not change anything bellow this line unless you know what you are doing
#-----------------------------------------------------------------------------

#
# Main
#     
from protocol_kerdensom import *

if __name__ == '__main__':
   protocolMain(ProtKerdensom)
