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
# {file}(*.xmd){validate}(PathExists) Input metadata:
InMetadata = ''

# Number of parts
"""The metadata will be splitted into this number of parts"""
Nparts=2

# {list_combo}(Do not sort, image name, micrograph name) Sort output by:
SortBy='Do not sort'

# {expert}Remove disabled images:
RemoveDisabled=True

# {expert}Split randomly:
RandomSplit=True

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_metadata_split import *
#        
# Main
#     
 
if __name__ == '__main__':
    # create preprocess_particles_class object
    protocolMain(ProtMetadataSplit)
