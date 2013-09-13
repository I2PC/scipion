#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of volumes 
# Author: Carlos Oscar, August 2013
#
# {begin_of_header}

# {eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input
#-----------------------------------------------------------------------------
# {file}(*.vol *.stk *.xmd){validate}(PathExists) Volumes to process:
""" Can be a density volume or a PDB file"""
InModel = ''

# Input voxel size (A/voxel):
InitialTs=1.0

#-----------------------------------------------------------------------------
# {section} Preprocess steps
#-----------------------------------------------------------------------------
# Change hand
""" Change hand by applying a mirror along X """
DoChangehand=False

# Randomize phases
""" Randomize phases beyond a certain frequency """
DoRandomize=False

# {condition}(DoRandomize==True) Maximum Resolution
""" Angstroms """
MaxResolutionRandomize=40

# Low pass filter
""" Low pass filter the volume """
DoFilter=False

# {condition}(DoFilter==True) Maximum Resolution
""" Angstroms """
MaxResolution=40

# Symmetrize
""" Symmetrize the input model """
DoSymmetrize=False

# {condition}(DoSymmetrize==True) Symmetry group
""" See [http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry]
    for a description of the symmetry groups format
    If no symmetry is present, give c1
"""
Symmetry='c1'

# {condition}(DoSymmetrize==True){list_combo}(average,sum) Aggregation mode
""" Symmetrized volumes can be averaged or summed.
"""
SymmetryAggregation='average'

# Spherical Mask
""" Remove all pixels beyond a certain radius """
DoMask=False

# {condition}(DoMask==True){wizard}(wizardSetMaskRadiusPreprocess) Mask Radius
""" In pixels. Set to -1 for half of the size of the volume """
MaskRadius=-1

# Adjust gray values
""" Adjust input gray values so that it is compatible with a set of projections """
DoAdjust=False

# {condition}(DoAdjust==True) Set of projections
""" Metadata with a set of images to which the model should conform. The set of images should have the
final pixel size and the final size of the model """
SetOfImages=''

# Normalize background
""" Set background to have zero mean and standard deviation 1 """
DoNormalize=False

# {condition}(DoNormalize==True){wizard}(wizardSetMaskRadiusPreprocess2) Mask Radius
""" In pixels. Set to -1 for half of the size of the volume """
MaskRadiusNormalize=-1

# Threshold
""" Remove voxels below a certain value """
DoThreshold=False

# {condition}(DoThreshold==True) Threshold value
""" Grey value below which all voxels should be set to 0 """
Threshold=0

# Segment
""" Separate the molecule from its background """
DoSegment=False

# {condition}(DoSegment==True) {list_combo}(Voxel mass, Aminoacid mass, Dalton mass, Automatic) Segmentation type
""" Type of segmentation """
SegmentationType='Automatic'

# {condition}(DoSegment==True and SegmentationType!="Automatic") Molecule Mass
""" In automatic segmentation, set it to -1"""
SegmentationMass=-1

#-----------------------------------------------------------------------------
# {section} Output
#-----------------------------------------------------------------------------
# Final voxel size (A/voxel):
FinalTs=1.0

# Final box size (voxels):
""" Set to -1 for automatic estimation  """
FinalSize=-1

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_preprocess_volume import *
#        
# Main
#     
 
if __name__ == '__main__':
    # create preprocess_particles_class object
    protocolMain(ProtPreprocessVolumes)
