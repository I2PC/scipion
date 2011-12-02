#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for random conical tilt reconstruction
#
# Example use:
# ./xmipp_protocol_rct.py
#
#  Author: Sjors Scheres, March 2007
# Updated: J.M. de la Rosa Trevin October 2011
#
# {begin_of_header}

#{eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input
#-----------------------------------------------------------------------------

# {run}(extract_particles) Particles extraction RUN
""" 
Select a previous run of Particles extraction protocol.
In the run directory should be a file named <tilted_pairs.xmd>
containing the pairs of untilted and tilted images 
(labels MDL_IMAGE and MDL_IMAGE_TILTED respectively)
"""
ExtractionRun = ''

# {file}(tilted_pairs.xmd){validate}(PathExists) Input tilt pairs metadata:
""" 
Provide a metadata containing the pairs of images that make up your data set.
This metadata should contain at least the labels MDL_IMAGE and MDL_IMAGE_TILTED
The filenames should be relative to the <ProjectDir> where you are running the <Protocols>
If you have images outside the <ProjectDir> you should import them first.
"""
TiltPairsMd = "tilted_pairs.xmd"

# {file}(result*classes.xmd){validate}(PathExists) 2D Classification metadata:
"""
You should provide a metadata where all images in the dataset are
grouped into 2D classes. For more details about this file, see:
[http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/ProtocolsIntercomunication]
"""
ClassifMd = ''

# Which classes do you want to reconstruct?
"""
Select which 2D classes do you want to reconstruct from
the previous classification run. 
Use comma-separated lists, with the "-" sign to indicate ranges
Example: <1,3,6-9,12>
For each class that you select a volume will be reconstructed
with the information of alignment of each untilted image and
the corresponding tilted pair.
"""
SelectedClasses = ''

#------------------------------------------------------------------------------------------------
# {section}{has_question} Prepare images
#------------------------------------------------------------------------------------------------

# Prepare local copies of all images?
""" 
This will make local copies of all untilted images and generate corresponding selfiles. 
This has to be done at least once.
"""
DoImagePreparation=True
# Set untilted image headers?
"""
This will re-align the untilted particles and set the RCT angles correctly
"""
DoUntiltedHeaders=True

# Set tilted image headers?
""" This will center the tilted particles and set the RCT angles correctly
"""
DoTiltedHeaders=True

# Maximum allowed shift for tilted particles (pixels):
""" 
Particles that shift more will be discarded. A value larger than the 
image size will not discard any particle.
"""
CenterMaxShift=999

# {expert} Additional alignment parameters
"""
This are additional parameters for the program: <xmipp_image_align_tilt_pairs>
For example:
<--skip_stretching> will skip the cosine-stretching prior to centering
<--skip_centering>  will skip the entire centering, so that only the RCT 
                    angles will be set.
<--force_x_zero>    will force the shift in the X direction to be zero, 
                    and will only center in the Y direction
For more details see: 
[http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Image_align_tilt_pairs_v3]
"""
AlignTiltPairsAdditionalParams = ""

#------------------------------------------------------------------------------------------------
# {section}{has_question} Reconstruction
#------------------------------------------------------------------------------------------------

# Reconstruct 3D-classes
DoReconstruct = True

# {list}(fourier,art,wbp) Reconstruction method
""" 
The <art> method is very slow, and the relaxation parameter (-l ) needs 
to be carefully tuned. The <wbp> and <fourier> methods are much faster, 
but <wbp> may give significantly worse results. 
Therefore, the <fourier> method is the recommended one.
"""
ReconstructMethod = 'fourier'

# {expert} Additional reconstruction parameters
"""
For <fourier> see: [http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Fourier]
For <art> see: [http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Art]
For <wbp> see: [http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Wbp]
"""
ReconstructAdditionalParams = ""

# Filter reconstructed volumes?
""" 
Filtering may be useful to remove noise, especially when few particles 
contribute to the reconstruction.
"""
DoLowPassFilter = True

# {condition}(DoLowPassFilter) Resolution of the low-pass filter (Ang):
LowPassFilter = 50

# {condition}(DoLowPassFilter) Pixel size (Ang):
PixelSize = 5.6

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------

from protocol_rct import *
#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtRCT)
