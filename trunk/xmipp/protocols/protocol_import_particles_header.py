#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based importing single-particles: 
#
# Author: Carlos Oscar, August 2011
#
# {begin_of_header}

# {eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Importing parameters
#-----------------------------------------------------------------------------
# {file}(*.xmd *.stk *.sel *.ctfdat){validate}(PathExists) Stack or Selfile
InputFile=''

# {validate}(NonEmpty) Family name
""" Name of the family with which they will be imported """
Family="Imported"

# {validate}(IsFloat) Sampling rate (A/pixel)
SamplingRate = 1

# Copy images
""" Copy images from their original location """
DoCopy=False

# Import all images in input
ImportAll=True

#-----------------------------------------------------------------------------
# {section}{condition}(ImportAll) Preprocessing parameters
#-----------------------------------------------------------------------------
# 1) Invert contrast
""" Invert the contrast if your particles are black over a white background. """
DoInvert=False

# 2) Dust particles removal
""" Sets pixels with unusually large values to random values from a Gaussian with zero-mean and unity-standard deviation.
    It requires a previous normalization, i.e., Normalization must be set to Yes.
"""
DoRemoveDust=False

# {expert} Threshold for dust removal:
""" Pixels with a signal higher or lower than this value times the standard deviation of the image will be affected. For cryo, 3.5 is a good value.
    For high-contrast negative stain, the signal itself may be affected so that a higher value may be preferable.
"""
DustRemovalThreshold=3.5

# 3) Normalize (Recommended)
""" It subtract a ramp in the gray values and normalizes so that in the background
    there is 0 mean and standard deviation 1 """
DoNorm=True

# {expert} Background radius
"""Pixels outside this circle are assumed to be noise and their stddev is set to 1.
   Radius for background circle definition (in pix.).
   If this value is 0, then the same as the particle radius is used. """
BackGroundRadius=0

#-----------------------------------------------------------------------------
# {section}{condition}(not ImportAll) Subset parameters
#-----------------------------------------------------------------------------
# {list_combo}(Random subset,First particles) Subset mode
SubsetMode="First particles"

# Import the first N particles
""" A subset with the first N particles will be generated """
Nsubset=500

# {eval} expandParallel(threads=0, mpi=8, hours=12, condition="ImportAll")

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_import_particles import *
#        
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtImportParticles)