#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of x-ray tomograms: 
#
# For each tomogram given, this script will perform 
# the requested operations below.
# For each tomogram a subdirectory will be created
#
# Author: Joaquin Oton, October 2013
#
#------------------------------------------------------------------------------------------------
# {begin_of_header}

#{eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input tomograms
#-----------------------------------------------------------------------------
# {run}(xray_fast_align) Align Tomograms Run
ImportRun = ''

# {list_combo}(imod, tomo3d) Reconstruction program
"""Tomo series from Mistral microscope (Alba synchrotron) are expected to be
 in nexus format in hdf5 files, while from U41-TXM line (Bessy) you only set 
 the data folder and the indexes for initial and final images for both tomogram 
 and flatfields."""
recProgram = "imod"

# {validate}(IsInt) Thickness
thickness = ''

# {expert} Contrast inversion
"""Before reconstruct invert the contrast of the images. It is mandatory to 
   recover real absorption coefficient volumes."""
DoInvertContrast = True

# {eval} expandParallel(threads=0,hours=6,mpi=1)

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_xray_reconstruct import *
if __name__ == '__main__':
    protocolMain(ProtXrayRecon)
