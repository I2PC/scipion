#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of volumes 
# Author: Carlos Oscar, August 2013
#
# {begin_of_header}

# {eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Create mask
#-----------------------------------------------------------------------------
# {list_combo}(Volume, Geometry, Binary mask file) Mask Source
""" Choose volume for thresholding or segmenting a volume. For geometrical masks
    choose geometry. For a binary mask already generated use binary mask file"""
MaskSource="Volume"

# {file}(*.vol){validate}(PathExists){condition}(MaskSource=="Volume") Input volume:
""" Density volume"""
InModel = ''

# {list_combo}(Threshold, Segment){condition}(MaskSource=="Volume") Operation
VolumeOperation="Threshold"

# {condition}(MaskSource=="Volume" and VolumeOperation=="Threshold") Threshold value
""" Grey value below which all voxels should be set to 0 """
Threshold=0

# {condition}(MaskSource=="Volume" and VolumeOperation=="Segment") {list_combo}(Voxel mass, Aminoacid mass, Dalton mass, Automatic) Segmentation type
""" Type of segmentation """
SegmentationType='Automatic'

# {condition}(MaskSource=="Volume" and VolumeOperation=="Segment" and SegmentationType!="Automatic") Molecule Mass
""" In automatic segmentation, set it to -1"""
SegmentationMass=-1

# {condition}(MaskSource=="Volume" and VolumeOperation=="Segment" and SegmentationType!="Automatic") Voxel size (A/vox)
Ts=1

# {condition}(MaskSource=="Geometry") Mask size
""" In voxels """
MaskSize=128

# {list_combo}(Sphere, Box, Crown, Cylinder, Gaussian, Raised cosine, Raised crown){condition}(MaskSource=="Geometry") Mask type
SimpleOperation="Sphere"

# {condition}(MaskSource=="Geometry" and [SimpleOperation=="Sphere" or SimpleOperation=="Cylinder"]) Radius
""" In pixels. """
Radius=50

# {condition}(MaskSource=="Geometry" and SimpleOperation=="Box") Box size
""" In pixels. """
BoxSize=50

# {condition}(MaskSource=="Geometry" and [SimpleOperation=="Crown" or SimpleOperation=="Raised cosine" or SimpleOperation=="Raised crown"]) Inner radius
""" In pixels. """
InnerRadius=40

# {condition}(MaskSource=="Geometry" and [SimpleOperation=="Crown" or SimpleOperation=="Raised cosine" or SimpleOperation=="Raised crown"]) Outer radius
""" In pixels. """
OuterRadius=50

# {condition}(MaskSource=="Geometry" and SimpleOperation=="Cylinder") Height
""" In pixels. """
Height=50

# {condition}(MaskSource=="Geometry" and SimpleOperation=="Gaussian") Sigma
""" In pixels """
Sigma=20

# {condition}(MaskSource=="Geometry" and SimpleOperation=="Raised crown") Pixel width
""" In pixels """
PixelWidth=5

# {condition}(MaskSource=="Binary mask file"){file}(*.vol){validate}(PathExists) Binary mask
""" File """
BinaryMask=""

#-----------------------------------------------------------------------------
# {section} Postprocess mask
#-----------------------------------------------------------------------------
# Symmetrize
""" Symmetrize mask """
DoSymmetrize=False

# {condition}(DoSymmetrize==True) Symmetry group
""" See [http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry]
    for a description of the symmetry groups format
    If no symmetry is present, give c1
"""
Symmetry='c1'

# Apply morphological operation
DoMorphological=False

# {list_combo}(dilation, erosion, closing, opening){condition}(DoMorphological) Operation
MorphologicalOperation='dilation'

# {expert}{condition}(DoMorphological) Structural element size
ElementSize=1

# Invert mask
DoInvert=False

# Smooth mask
DoSmooth=False

# {condition}(DoSmooth) Gaussian sigma
""" Smoothing is performed by convolving the mask with a Gaussian of this sigma (in pixels)"""
SigmaConvolution=2

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_create_volume_mask import *
#        
# Main
#     
 
if __name__ == '__main__':
    # create preprocess_particles_class object
    protocolMain(ProtCreateVolumeMask)
