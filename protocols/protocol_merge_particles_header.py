#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based merging of different particle sets 
#
# Author: Roberto Marabini        (roberto@cnb.csic.es)     July 2013
#         J. M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
#

# {begin_of_header}

# {eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input images
#-----------------------------------------------------------------------------

# {file}(images*.xmd){validate}(PathExists) First set of images:
""" Select a set of particle images to merge with another set. """
InputImages1 = ''

# {file}(images*.xmd){validate}(PathExists) Second set of images:
""" Select a second set to be merged. """
InputImages2 = ''

# {hidden} Protocol help
Usage = """
Merge two particle sets (usually images.xmd files) to produce
another set (stored in images.xmd on the protocol working dir)
"""

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_merge_particles import *
#        
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtMergeParticles)