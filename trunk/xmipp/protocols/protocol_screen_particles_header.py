#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based screening single-particles: 
# Author: Carlos Oscar, October 2011
#
# {begin_of_header}

#{include} inc_comment_run.py

#-----------------------------------------------------------------------------
# {section} Extracting parameters
#-----------------------------------------------------------------------------
# {file}{validate}(PathExists) Input images
""" Stack file or metadata file with a selection file"""
InputFile=''

#------------------------------------------------------------------------------------------
# {section} Parallelization
#------------------------------------------------------------------------------------------
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
from protocol_screen_particles import *
#        
# Main
#     
 
if __name__ == '__main__':
       # create preprocess_particles_class object
    protocolMain(ProtScreenParticles)