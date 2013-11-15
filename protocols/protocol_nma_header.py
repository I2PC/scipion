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
PseudoAtomRadius=1

# {expert} Volume approximation error (%):
"""This value is a percentage (between 0.001 and 100) specifying how fine you want to
approximate the EM volume by the pseudoatomic structure. Lower values imply lower approximation error, and consequently, more pseudoatoms."""
PseudoAtomTarget=5

#------------------------------------------------------------------------------------------------
# {section} Normal Mode Analysis
#------------------------------------------------------------------------------------------------
# Number of modes:
"""The maximum number of modes allowed by the method for atomic normal mode analysis is 6 times the number of RTB blocks and for pseudoatomic normal mode analysis 3 times the number of pseudoatoms. However, the protocol allows only up to 200 modes as 20-100 modes are usually enough. The number of modes given here should be below the minimum between these two numbers."""
NumberOfModes=20

# {list}(Absolute,Relative) Cut-off mode:
CutoffMode="Relative"

# {condition}(CutoffMode=="Absolute") Cut-off distance (A):
"""Atoms or pseudoatoms beyond this distance will not interact"""
Rc=8

# {condition}(CutoffMode=="Relative") Cut-off percentage:
"""The interaction cutoff distance is calculated as the distance below which is this percentage of interatomic or interpseudoatomic distances. Atoms or pseudoatoms beyond this distance will not interact. """
RcPercentage=95

# {expert}{condition}(StructureType=="PDB") Number of residues per RTB block:
"""This is the RTB block size for the RTB NMA method. 
When calculating the normal modes, aminoacids are grouped into blocks of this size that are moved translationally and rotationally together """
RTBblockSize=10

# {expert}{condition}(StructureType=="PDB") Interaction force constant
"""This is the RTB interaction force constant for the RTB NMA method. 
If it increases, then the structure will be more rigid. """
RTBForceConstant=10.0

# {expert} Threshold on collectivity
""" Collectivity degree is related to the number of atoms or pseudoatoms that are affected by the mode, and it is normalized between 0 and 1. Modes below this threshold are deselected in the modes metadata file. Set to 0 for no deselection. You can always modify the selection manually after the modes metadata file is created. The modes metadata file can be used with Flexible fitting protocol. Modes 1-6 are always deselected as they are related to rigid-body movements.  """
CollectivityThreshold=0.15

#------------------------------------------------------------------------------------------------
# {section} Animation
#------------------------------------------------------------------------------------------------
# Amplitude:
Amplitude=50

#{expert} Number of frames:
NFrames=10

#{expert}{condition}(StructureType=="EM") Downsample pseudoatomic structure:
""" Downsample factor 2 means removing one half of the pseudoatoms """
Downsample=1

#{expert}{condition}(StructureType=="EM") Pseudoatom mass threshold:
""" Remove pseudoatoms whose mass is below this threshold. This value should be between 0 and 1"""
PseudoAtomThreshold=0

#------------------------------------------------------------------------------------------------
# {section}{visualize} Visualization
#------------------------------------------------------------------------------------------------

# {view}{condition}(StructureType=="EM") Display pseudoatom representation
DisplayPseudoAtom=True

# {view}{condition}(StructureType=="EM") Display pseudoatom approximation
DisplayPseudoApproximation=True

# {view} Open list of modes
DisplayModes=True

# {view} Plot max distance profile
DisplayMaxDistanceProfile=True

# {view} Plot distance profile
DisplayDistanceProfile=True

# Open specific mode
DisplaySingleMode="7"

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#		
# Main
#     
from protocol_nma import *

if __name__ == '__main__':
    protocolMain(ProtNMA)
