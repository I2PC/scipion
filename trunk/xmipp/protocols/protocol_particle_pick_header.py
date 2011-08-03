#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based manual particle picking
#
# For each micrograph in the PreprocessingDir, this program will launch
# the xmipp_mark program 
# A graphical interface exists to identify micrographs that have been finished
#
# Author: Sjors Scheres, March 2007
# Author: Carlos Oscar Sorzano, June, 2011
#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
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

#------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------
# Run name
""" Working directory for this protocol 
"""
RunName = "particles_001"

# {run}(preprocess_micrographs) Preprocessing Micrographs run
""" Directory with the preprocessing (output of the Preprocessing Micrographs protocol)
"""
PreprocessingRun = "Preprocessing/micrographs_001"

# Perform automatic particle picking
""" Perform automatic particle picking """
AutomaticPicking = True

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
#------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
#------------------------------------------------------------------------------------------------

from protocol_particle_pick import *
#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtParticlePicking)
