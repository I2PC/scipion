#!/usr/bin/env python
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
Behavior = "Restart"

#------------------------------------------------------------------------------------------------
# {section} Common line parameters
#------------------------------------------------------------------------------------------------
# {file}{validate}(PathExists) Selfile with the input class averages:
""" This selfile points to the stack or metadata containing your images 
"""
InSelFile=''

# Particle radius
Radius=25
""" Maximum radius of the particle
"""
# Symmetry
""" c1=No symmetry; c2, c3, ...
"""
Symmetry='c1'

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
QueueHours = 96

# {hidden} Show expert options
"""If True, expert options will be displayed
"""
ShowExpertOptions = False

# Number of MPI processes to use
NumberOfMpi = 3

# {hidden} Show expert options
"""If True, expert options will be displayed
"""
ShowExpertOptions = False
#-----------------------------------------------------------------------------
# {end_of_header} do not change anything bellow this line unless you know what you are doing
#-----------------------------------------------------------------------------

from protocol_commonlines import *
  
if __name__ == '__main__':
   protocolMain(ProtCommonLines)
