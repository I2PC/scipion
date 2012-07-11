#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based extracting single-particles: 
#  - phase flipping
#  - extraction of particles
#  - normalization
#  - sort_junk
#
# It is assumed that you have already ran the preprocess_micrographs protocol,
#  and that you have picked the particles for each micrograph
#
# Author: Carlos Oscar, August 2011
#
# {begin_of_header}

# {eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Extracting parameters
#-----------------------------------------------------------------------------
# {list} (original, same as picking, other) Select downsampling factor:
DownsampleType = 'same as picking'

# {condition}(DownsampleType == 'other') Extraction downsampling factor
""" 
This factor is always referred to the original sampling rate.
You may use independent downsampling factors for extracting the 
particles, picking them and estimating the CTF. All downsampling 
factors are always referred to the original sampling rate, and 
the differences are correctly handled by Xmipp.
"""
DownsampleFactor = 2

# {run}(particle_pick,particle_pick_supervised, particle_pick_auto) Particle picking run
PickingRun = ''

# {wizard}(wizardChooseFamilyToExtract) Family 
Family = ''

# Particle box size
""" In pixels. The box size is the size of the boxed particles,
actual particles may be smaller than this. """
ParticleSize = 0

# Phase flipping (Recommended)
""" Use the information from the CTF to compensate for phase reversals."""
DoFlip = True

# Invert contrast
""" Invert the contrast if your particles are black over a white background. """
DoInvert = False

# Normalize (Recommended)
""" 
It subtract a ramp in the gray values and normalizes so that in the 
background there is 0 mean and standard deviation 1 """
DoNorm = True

# {expert}{condition}(DoNorm) Background radius
"""
Pixels outside this circle are assumed to be noise and their stddev 
is set to 1. Radius for background circle definition (in pix.).
If this value is 0, then half the box size is used. """
BackGroundRadius = 0

# {condition}(DoNorm) Dust particles removal
""" 
Sets pixels with unusually large values to random values from a Gaussian
with zero-mean and unity-standard deviation. 
"""
DoRemoveDust = False

# {expert}{condition}(DoRemoveDust) Threshold for dust removal:
""" 
Pixels with a signal higher or lower than this value times the standard 
deviation of the image will be affected. For cryo, 3.5 is a good value.
For high-contrast negative stain, the signal itself may be affected so 
that a higher value may be preferable.
"""
DustRemovalThreshold = 3.5



# {eval} expandParallel(threads=0, hours=12, mpi=8)
#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_extract_particles import *
#        
# Main
#     
 
if __name__ == '__main__':
    # create preprocess_particles_class object
    protocolMain(ProtExtractParticles)