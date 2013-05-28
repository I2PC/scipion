#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Normal Mode analysis of atomic and EM structures
#
# Author: Carlos Oscar Sanchez Sorzano, May 2013
#         Slavica Jonic
#

# {begin_of_header}

# {eval} expandCommentRun()

#------------------------------------------------------------------------------------------------
# {section} Volume to analyze
#------------------------------------------------------------------------------------------------
# {file}{validate}(PathExists) Input structure:
InputStructure=''

# {list}(EM, PDB) Kind of structure
""" Choose PDB for atomic structures and EM for Electron Microscopy Volumes 
"""
StructureType='EM'

# {expert}{condition}(StructureType=="EM")Volume approximation error
"""This value is a percentage (between 0 and 1) specifying how fine you want to
approximate the EM volume by the pseudoatomic structure. Lower values imply lower approximation error,
and consequently, more atoms."""
PseudoAtomTarget=0.02

# {eval} expandParallel(threads=0)

#------------------------------------------------------------------------------------------------
# {section}{visualize} Visualization
#------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#		
# Main
#     
from protocol_nma import *

if __name__ == '__main__':
    protocolMain(ProtNMA)
