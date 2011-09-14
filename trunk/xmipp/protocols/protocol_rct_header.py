#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for random conical tilt reconstruction
#
# Example use:
# ./xmipp_protocol_rct.py
#
# Author: Sjors Scheres, March 2007
#
# {begin_of_header}

#------------------------------------------------------------------------------------------
# {section}{has_question} Comment
#------------------------------------------------------------------------------------------
# Display comment
DisplayComment = False

# {text} Write a comment:
""" 
Describe your run here...
"""
#-----------------------------------------------------------------------------
# {section} Global parameters
#-----------------------------------------------------------------------------
# Run name:
""" This will identify your protocol run. It need to be unique for each protocol. You could have run1, run2 for protocol X, but not two
run1 for it. This name together with the protocol output folder will determine the working dir for this run.
"""
RunName = "run_001"

# {file} Selfile with all untilted images:
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
UntiltedSelFile="all_images_untilted.sel"
# {file} Selfile with all tilted images:
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
TiltedSelFile="all_images_tilted.sel"

#------------------------------------------------------------------------------------------------
# {section} Previous ML2D classification (WITHOUT INCLUDING MIRRORS!)
#------------------------------------------------------------------------------------------------
# {dir} Classification directory
"""Directory of previous ML2D-classification on the untilted images"""
PreviousDirML2D="ML2D/ML12ref"
# {expert} Rootname for ML2D run
""" Only provide something if the ml2d output rootname was different from the default.
    This will never be the case when using the standardized protocols 
"""
PreviousML2DRoot=""
# Which of these classes do you want to reconstruct? 
""" Use comma-separated lists, with the "-" sign to indicate ranges
"""
SelectClasses="1,3,6-9,12"
#------------------------------------------------------------------------------------------------
# {section}{has_question}Prepare images
#------------------------------------------------------------------------------------------------
# Prepare local copies of all images?
""" This will make local copies of all untilted images and generate corresponding selfiles. This has to be done at least once.
"""
DoImagePreparation=True
# Set untilted image headers?
""" This will re-align the untilted particles and set the RCT angles correctly
"""
DoUntiltedHeaders=True
# Set tilted image headers?
""" This will center the tilted particles and set the RCT angles correctly
"""
DoTiltedHeaders=True
# Maximum allowed shift for tilted particles (in pixels):
""" Particles that shift more will be discarded. A value larger than the image size will not discard any particle
"""
CenterMaxShift=999
# {expert} Additional parameters for the align_tilt_pairs program
"""  For example:
    -skip_stretching will skip the cosine-stretching prior to centering
    -skip_centering  will skip the entire centering, so that only the RCT angles will be set.
    -force_x_zero    will force the shift in the X direction to be zero, and will only center in the Y direction
"""
AlignTiltPairsAdditionalParams=""
#------------------------------------------------------------------------------------------------
# {section}{has_question} Reconstruction
#------------------------------------------------------------------------------------------------
# Reconstruct 3D-classes
DoReconstruct=True
# {list}(fourier,art,wbp) Reconstruction method
""" The art method is very slow, and the relaxation parameter (-l ) needs to be carefully tuned. 
    The wbp and fourier methods are much faster, but wbp may give significantly worse results.
    Therefore, the fourier method is the recommended one.
"""
ReconstructMethod='fourier'
# {expert} Additional reconstruction parameters
""" For fourier see: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Fourier
    For art see: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Art
    For wbp see: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Wbp
"""
ReconstructAdditionalParams=""
#------------------------------------------------------------------------------------------------
# {section}{has_question} Low-pass filter 
#------------------------------------------------------------------------------------------------
# Filter reconstructed volumes
""" Filtering may be useful to remove noise, especially when few particles contribute to the reconstruction.
"""
DoLowPassFilter=True
# Resolution of the low-pass filter (in Angstroms):
LowPassFilter=50
# Pixel size (in Angstroms):
PixelSize=5.6

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------

from protocol_rct import *
#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtRCT)
