#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of volumes 
# Author: Carlos Oscar, August 2013
#
# {begin_of_header}

# {eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input
#-----------------------------------------------------------------------------
# {file}(*.pdb){validate}(PathExists) Input model:
""" PDB file"""
InModel = ''

# Final voxel size (A/voxel):
FinalTs=1.0

# Final box size (voxels):
""" Set to -1 for automatic estimation  """
FinalSize=-1

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_convert_pdb import *
#        
# Main
#     
 
if __name__ == '__main__':
    # create preprocess_particles_class object
    protocolMain(ProtConvertPDB)
