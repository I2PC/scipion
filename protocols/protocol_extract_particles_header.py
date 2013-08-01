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
# {section} Input
#-----------------------------------------------------------------------------
# {list_combo} (original, same as picking, other) Downsampling type:
DownsampleType = 'same as picking'

# {condition}(DownsampleType == 'other') Downsampling factor:
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

# {wizard}(wizardChooseSizeToExtract) Particle box size
""" In pixels. The box size is the size of the boxed particles,
actual particles may be smaller than this. """
ParticleSize = 0

# {expert}{list_combo}(None, MaxZscore, Percentage) Automatic particle rejection
""" How to automatically reject particles. It can be none (no rejection), maxZscore (reject a particle
if its Zscore is larger than this value), Percentage (reject a given percentage in each one of the screening criteria). """
RejectionMethod='None'

# {expert}{condition}(RejectionMethod=="MaxZscore") Maximum Zscore
MaxZscore=3

# {expert}{condition}(RejectionMethod=="Percentage") Percentage (%)
Percentage=5

# {eval} expandParticlesPreprocess(allowFlip=True)

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