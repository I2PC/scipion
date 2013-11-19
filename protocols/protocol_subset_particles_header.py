#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# Author: Carlos Oscar, October 2011
#
# {begin_of_header}

#{include} inc_comment_run.py

# {eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Subset parameters
#-----------------------------------------------------------------------------
#{hidden} Usage of the program
Usage = """Generate a metadata that is a subset of the full image set using the itemIds of the subset
"""

# {file}(images*.xmd){validate}(PathExists) Full image set
"""This is the complete set of images. It must be a metadata with an itemId column. Those items whose id are also in the
subset metadata will be selected."""
InputFile=''

# {file}(images*.xmd){validate}(PathExists) Subset
"""This metadata defines which items will be selected from the full image set. Those items in the full set whose items
are in the list of itemsId given by this file will be selected"""
SubsetFile=''

#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_subset_particles import *
#        
# Main
#     
 
if __name__ == '__main__':
       # create preprocess_particles_class object
    protocolMain(ProtSubsetParticles)
