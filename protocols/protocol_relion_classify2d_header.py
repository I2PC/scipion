#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Relion-based 3D classification (relion 1.2)
#
# Author: J. M. de la Rosa Trevin     (jmdelarosa@cnb.csic.es) Feb 2014
#
#------------------------------------------------------------------------------------------------
# {begin_of_header}

# {eval} expandCommentRun()

# {eval} expandRelion(classify=True, is2D=True)

# {eval} expandParallel()

#------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
#------------------------------------------------------------------------------------------------

from protocol_relion_classify2d import *

if __name__ == '__main__':
    protocolMain(ProtRelion2DClassifier)

