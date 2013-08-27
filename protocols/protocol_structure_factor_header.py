#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# Author: Carlos Oscar, August 2013
#
# {begin_of_header}

# {eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input volume
#-----------------------------------------------------------------------------
# {file}(*.vol){validate}(PathExists) Input volume:
""" Density volume"""
InModel = ''

# Voxel size (A/vox)
Ts=1

#------------------------------------------------------------------------------------------------
# {section}{visualize} Visualization
#------------------------------------------------------------------------------------------------

# Display Structure factor
DisplayStructureFactor=True

# Display Guinier plot
DisplayGuinier=True

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_structure_factor import *
#        
# Main
#     
 
if __name__ == '__main__':
    # create preprocess_particles_class object
    protocolMain(ProtStructureFactor)
