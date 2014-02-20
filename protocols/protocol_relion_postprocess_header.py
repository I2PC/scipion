#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol header paramas for calling relion_postprocess program
#
# Author:  J. M. de la Rosa Trevin (jmdelarosa@cnb.csic.es) Sept 2011
#
#------------------------------------------------------------------------------------------------
# {begin_of_header}

# {eval} expandCommentRun()
# {cite}
CiteML3D = """
for ML3D:  Scheres et al. (2007) Nature Methods, 4, 27-29
for MLF3D: Scheres et al. (2007) Structure, 15, 1167-1177
"""

#------------------------------------------------------------------------------------------
# {section} Input
#------------------------------------------------------------------------------------------
# {file}(relion_model.star){validate}(PathExists) Input refine file:
""" 
Select the <relion_model.star> file of the 3D auto-refine Relion run
that you want to used for proceed with the post-process.
"""
ModelStar = ""

#------------------------------------------------------------------------------------------
# {section} Masking
#------------------------------------------------------------------------------------------
# Perform automated masking?
"""
Perform automated masking, based on a density threshold
Relion param: <--auto_mask>
"""
AutoMask = False

# {condition}(not AutoMask) Initial mask threshold
"""
Density at which to threshold the map for the initial seed mask
Relion param: <--inimask_threshold>
"""
InimaskThreshold = 0.02

# {condition}(not AutoMask) Mask pixels extension (pix)
"""
Number of pixels to extend the initial seed mask
Relion param: <--extend_inimask>
"""
ExtendInimask = 3.

# {condition}(not AutoMask) {condition}(not AutoMask) Mask edge width (pix)
"""
Width for the raised cosine soft mask edge (in pixels)
Relion param: <--width_mask_edge>
"""
WidthMaskEdge = 6.

# {condition}(not AutoMask) Mask file
"""
Filename of a user-provided mask 
(1=protein, 0=solvent, all values in range [0,1])
Relion param: <--mask>
"""
Mask = ""

#------------------------------------------------------------------------------------------
# {section} Sharpening
#------------------------------------------------------------------------------------------

# {file}(*.star) MTF-curve file:
"""
User-provided STAR-file with the MTF-curve of the detector
Relion param: <--mtf>
"""
Mtf = ""

# Determine B-factor automatically?
"""
Perform automated B-factor determination.
<(Rosenthal and Henderson, 2003)>
Relion param: <--auto_bfac>
"""
AutoBfac = False

# {condition}(AutoBfac) B-factor lowest resolution (A)
"""
Lowest resolution (in A) to include in fitting of the B-factor
Relion param: <--autob_lowres>
"""
AutobLowres = 10.

# {condition}(AutoBfac) B-factor highest resolution (A)
"""
Highest resolution (in A) to include in fitting of the B-factor
Relion param: <--autob_highres>
"""
AutobHighres = 0.

# {condition}(not AutoBfac) Provide B-factor
"""
User-provided B-factor (in A^2) for map sharpening, e.g. -400
Relion param: <--adhoc_bfac>
"""
AdhocBfac = 0.

#------------------------------------------------------------------------------------------
# {section} Filtering
#------------------------------------------------------------------------------------------
# Use FSC-weighting for sharpening?
"""
If set to <No>, do not use FSC-weighting in the sharpening process
<(Rosenthal and Henderson, 2003)>
Relion param: <--skip_fsc_weighting>
"""
UseFscWeighting = True

# Low-pass filter (A)
"""
Resolution (in Angstroms) at which to low-pass filter the final map 
(by default at final resolution)
Relion param: <--low_pass>
"""
LowPass = 0.

# {expert} Low-pass filter edge width
"""
Width of the raised cosine on the low-pass filter edge 
(in resolution shells)
Relion param: <--filter_edge_width>
"""
FilterEdgeWidth = 2

# {expert} Randomize phases threshold
"""
Randomize phases from the resolution where FSC drops below this value
Relion param: <--randomize_at_fsc>
"""
RandomizeAtFsc = 0.8

# {expert} Verbosity
"""
Output verbosity.
Relion param: <--verb>
"""
Verb = 1

#------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
#------------------------------------------------------------------------------------------------

from protocol_relion_postprocess import *

if __name__ == '__main__':
    protocolMain(ProtRelionPostProcess)

