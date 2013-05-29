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

# {list}(EM, PDB) Kind of structure:
""" Choose PDB for atomic structures and EM for Electron Microscopy Volumes 
"""
StructureType='EM'

#------------------------------------------------------------------------------------------------
# {section}{condition}(StructureType=="EM") Pseudoatom conversion
#------------------------------------------------------------------------------------------------
# Sampling rate (A/voxel):
Sampling=1

# {list_combo}(No mask, Threshold, Binary mask) Mask mode:
MaskMode='No mask'

# {condition}(MaskMode=="Threshold") Threshold value:
Threshold=0.01

# {condition}(MaskMode=="Binary mask") Mask file:
MaskFile=""

#{expert} Pseudoatom radius (voxels):
""" Pseudoatoms are defined as Gaussians whose standard deviation is this value """
PseudoAtomRadius=2

# {expert} Volume approximation error (%):
"""This value is a percentage (between 0.001 and 100) specifying how fine you want to
approximate the EM volume by the pseudoatomic structure. Lower values imply lower approximation error,
and consequently, more atoms."""
PseudoAtomTarget=5

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
