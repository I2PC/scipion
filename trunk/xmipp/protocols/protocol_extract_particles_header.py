#!/usr/bin/env python
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
# {section} Run parameters
#-----------------------------------------------------------------------------
# Run name:
""" This will identify your protocol run. It need to be unique for each protocol. You could have run1, run2 for protocol X, but not two
run1 for it. This name together with the protocol output folder will determine the working dir for this run.
"""
RunName = "run_001"

# {list}(Resume, Restart) Run behavior
""" Resume from the last step, restart the whole process
"""
Behavior = "Resume"

#-----------------------------------------------------------------------------
# {section} Extracting parameters
#-----------------------------------------------------------------------------
# {run}(particle_pick,particle_pick_auto) Particle picking run
PickingRun=''

# {wizard}(wizardChooseFamilyToExtract) Family 
Family=''

# Particle size
""" In pixels """
ParticleSize=0

# Take Logarithm
""" Depending on your acquisition system you may have to take the logarithm
    or not in order to have a linear relationship between the gray values
    in the image and those in the volume """
DoLog=False 

# Phase flipping (Recommended)
""" Use the information from the CTF to compensate for phase reversals."""
DoFlip=True

# Invert contrast
""" Invert the contrast if your particles are black over a white background. """
DoInvert=False

# Normalize (Recommended)
""" It subtract a ramp in the gray values and normalizes so that in the background
    there is 0 mean and standard deviation 1 """
DoNorm=True

# {expert} Background radius
"""Pixels outside this circle are assumed to be noise and their stddev is set to 1.
   Radius for background circle definition (in pix.).
   If this value is 0, then the same as the particle radius is used. """
BackGroundRadius=0

# Dust particles removal
""" Sets pixels with unusually large values to random values from a Gaussian with zero-mean and unity-standard deviation.
    It requires a previous normalization, i.e., Normalization must be set to Yes.
"""
DoRemoveDust=False

# {expert} Threshold for dust removal:
""" Pixels with a signal higher or lower than this value times the standard deviation of the image will be affected. For cryo, 3.5 is a good value.
    For high-contrast negative stain, the signal itself may be affected so that a higher value may be preferable.
"""
DustRemovalThreshold=3.5

#------------------------------------------------------------------------------------------
# {section} Parallelization
#------------------------------------------------------------------------------------------
# Number of MPI processes
""" Set to 1 if you do not have MPI installed"""
NumberOfMpi = 8

# Submit to queue
"""Submit to queue
"""
SubmitToQueue = False

# {condition}(SubmitToQueue) Queue name
"""Name of the queue to submit the job
"""
QueueName = "default"

# {condition}(SubmitToQueue) Queue hours
"""This establish a maximum number of hours the job will
be running, after that time it will be killed by the
queue system
"""
QueueHours = 12

# {hidden} Show expert options
"""If True, expert options will be displayed
"""
ShowExpertOptions = False
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