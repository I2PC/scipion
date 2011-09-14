#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for building initial references by common lines
# This protocol is based on EMAN 1
#
# Example use:
# ./xmipp_protocol_commonlines.py
#
# Author:Carlos Oscar Sorzano, January 2011
#

# {begin_of_header}

#-----------------------------------------------------------------------------
# {section} Global parameters
#-----------------------------------------------------------------------------
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

# {file} Selfile or stack with the input images:
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
SelFileName='CL2D/classes8/class_level_00.stk'

#------------------------------------------------------------------------------------------------
# {section} Common lines parameters
#------------------------------------------------------------------------------------------------
# Particle radius
Radius=25
""" Maximum radius of the particle
"""
# Symmetry
""" c1=No symmetry; c2, c3, ...
"""
Symmetry='c1'

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
#-----------------------------------------------------------------------------
# {end_of_header} do not change anything bellow this line unless you know what you are doing
#-----------------------------------------------------------------------------

from protocol_commonlines import *

#
# Main
#     
if __name__ == '__main__':
   protocolMain(ProtCommonLines)
