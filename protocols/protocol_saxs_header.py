#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Author: Carlos Oscar Sanchez Sorzano, September 2013
#

# {begin_of_header}

# {eval} expandCommentRun()

#------------------------------------------------------------------------------------------------
# {section} Pseudoatom conversion
#------------------------------------------------------------------------------------------------
# {file}{validate}(PathExists) Input volume:
""" A SAXS curve will be computed from this volume. For doing so, it is first converted into
    pseudoatoms, and then we use crysol from atsas ([http://www.embl-hamburg.de/biosaxs/manuals/crysol.html])
    to calculate the theoretical SAXS curve that should be observed
"""
InputStructure=''

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
approximate the EM volume by the pseudoatomic structure. Lower values imply lower approximation error, and consequently, more pseudoatoms."""
PseudoAtomTarget=5

#------------------------------------------------------------------------------------------------
# {section} SAXS curve
#------------------------------------------------------------------------------------------------
# Experimental SAXS curve
""" This parameter is optional. If it is given the simulated SAXS curve will be compared to the
experimental one"""
ExperimentalSAXS=""

# {expert} Number of harmonics
""" See parameter lm [http://www.embl-hamburg.de/biosaxs/manuals/crysol.html] """ 
NumberOfHarmonics=18

# {expert} Maximum frequency
""" See parameter sm [http://www.embl-hamburg.de/biosaxs/manuals/crysol.html] """ 
MaxFreq=0.3

# {expert} Number of samples
""" See parameter ns [http://www.embl-hamburg.de/biosaxs/manuals/crysol.html] """ 
NumberOfSamples=256

# {expert} Other Crysol parameters
""" See [http://www.embl-hamburg.de/biosaxs/manuals/crysol.html] """ 
OtherCrysol=""


#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#		
# Main
#     
from protocol_saxs import *

if __name__ == '__main__':
    protocolMain(ProtSaxs)
