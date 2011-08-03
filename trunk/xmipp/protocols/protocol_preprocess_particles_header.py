#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of single-particles: 
#  - phase flipping
#  - extraction of particles
#  - normalization
#  - sort_junk
#
# It is assumed that you have already ran the preprocess_micrographs protocol,
#  and that you have picked the particles for each micrograph
#
# Example use:
# ./xmipp_preprocess_particles.py
#
# Author: Sjors Scheres, March 2007
#         Carlos Oscar, December 2010
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

#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Run name:
""" This will identify your protocol run. It need to be unique for each protocol. You could have run1, run2 for protocol X, but not two
run1 for it. This name together with the protocol output folder will determine the working dir for this run.
"""
RunName = "run_001"

# Delete working directory?
""" If TRUE the working directory will be deleted before run.
Set this option to TRUE if you want to start from scratch the same run
with previous parameters
"""
DoDeleteWorkingDir = False

# {run}(particle_pick) Particle picking run
PickingDir='ParticlePicking'

# {expert} Name for the output selfile:
""" This name should have extension .sel
"""
OutSelFile='all_images.sel'

#------------------------------------------------------------------------------------------------
# {section} Processing parameters
#------------------------------------------------------------------------------------------------
# Box size of the particles to extract (in pix.)
Size=80

# Do phase flipping?
DoFlip=True

# {expert} Take Logarithm?
DoLog=False 

# Invert contrast?
DoInvert=False

# {expert} Background radius
"""Pixels outside this circle are assumed to be noise and their stddev is set to 1.
   Radius for background circle definition (in pix.).
   If this value is 0, then the same as the particle radius is used. """
BackGroundRadius=0

# Perform dust particles removal?
""" Sets pixels with unusually large values to random values from a Gaussian with zero-mean and unity-standard deviation.
"""
DoRemoveDust=True

# {expert} Threshold for dust removal:
""" Pixels with a signal higher or lower than this value times the standard deviation of the image will be affected. For cryo, 3.5 is a good value. For high-contrast negative stain, the signal itself may be affected so that a higher value may be preferable.
"""
DustRemovalThreshold=3.5

#------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------
# Number of (shared-memory) threads?
""" This option provides shared-memory parallelization on multi-core machines.
It does not require any additional software, other than xmipp
"""
NumberOfThreads = 1

# Number of MPI processes to use
NumberOfMpi = 3

#------------------------------------------------------------------------------------------
# {section}{has_question} Queue
#------------------------------------------------------------------------------------------
# Submmit to queue
"""Submmit to queue
"""
SubmmitToQueue = True

# Queue name
"""Name of the queue to submit the job
"""
QueueName = "default"

# Queue hours
"""This establish a maximum number of hours the job will
be running, after that time it will be killed by the
queue system
"""
QueueHours = 72

# {hidden} Show expert options
"""If True, expert options will be displayed
"""
ShowExpertOptions = False
#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#  {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------

#        
# Main
#    
from protocol_preprocess_particles import *
 
if __name__ == '__main__':
       # create preprocess_particles_class object
    protocolMain(ProtPreprocessParticles)