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

# {run}(extract_particles,preprocess_particles) Particles extraction RUN
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

# {file}(classes*.xmd){validate}(PathExists) 2D Classification metadata:
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
# {section} Alignment parameters
#------------------------------------------------------------------------------------------------
# Thin Object
""" If the object is thin, then the tilted projections can be stretched to match the untilted projections"""
ThinObject=True

#{expert} Maximum allowed shift for tilted particles (pixels):
""" 
Particles that shift more will be discarded. A value larger than the 
image size will not discard any particle.
"""
CenterMaxShift=10

#{expert} Skip tilted translation alignment
""" If the tilted image quality is very low, then this alignment might result in poor estimates."""
SkipTiltedTranslations=False

#------------------------------------------------------------------------------------------------
# {section} Reconstruction
#------------------------------------------------------------------------------------------------
# {expert} Additional reconstruction parameters
"""
See: [http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Reconstruct_art_v31]
"""
ReconstructAdditionalParams = "-n 5 -l 0.01"

# Filter reconstructed volumes?
""" 
Filtering may be useful to remove noise, especially when few particles 
contribute to the reconstruction.
"""
DoLowPassFilter = True

# {condition}(DoLowPassFilter) Resolution of the low-pass filter (dig.freq):
LowPassFilter = 0.2

# {eval} expandParallel(threads=0,hours=12)

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------

from protocol_rct import *
#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtRCT)
