#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles, according to:
#
# Example use:
# ./xmipp_protocol_ml2d.py
#
#  Author:  Sjors Scheres, January 2008
# Updated:  J. M. de la Rosa Trevin July 2011
#
# {begin_of_header}
#{please_cite}
"""
for ML2D:  Scheres et al. (2005) J.Mol.Biol 348, 139-149
for MLF2D: Scheres et al. (2007) Structure 15, 1167-1177
"""
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
""" Resume from the last step, restart the whole process or continue at a given step or iteration
"""
Behavior = "Resume"

#-----------------------------------------------------------------------------
# {section} Input
#-----------------------------------------------------------------------------
# {file}{view}{validate}(PathExists)Metadata to be splitted
""" Metadata to be splitted
"""
InputMd = "all_images.xmd"

# Randomize metadata before split?
""" Perform randomization of items before split
"""
Randomize = True

# {condition}(Randomize)
# Random seed
""" Just testing conditions
"""
RandomSeed = 3

# {validate}(IsInt) Number of part
""" Number of new metadatas splitted from the input one
"""
NumberOfParts = 2

# {validate}(NonEmpty) Part prefix
Prefix = "part"

#------------------------------------------------------------------------------------------
# {section}{condition}(NumberOfParts>3)Parallelization 
#------------------------------------------------------------------------------------------
# Number of threads
""" This option provides shared-memory parallelization on multi-core machines.
It does not require any additional software, other than xmipp
"""
NumberOfThreads = 1

# Number of MPI processes
NumberOfMpi = 3

# Submit to queue ? 
"""Submit to queue
"""
SubmitToQueue = True

# {expert}{condition}(SubmitToQueue) Queue name
"""Name of the queue to submit the job
"""
QueueName = "default"

# {condition}(SubmitToQueue) {wizard}(wizardBrowse) Queue hours
"""This establish a maximum number of hours the job will
be running, after that time it will be killed by the
queue system
"""
QueueHours = 72

#------------------------------------------------------------------------------------------------
# {section}{visualize} Visualization
#------------------------------------------------------------------------------------------------
# Visualize the class averages of all iterations in matrix-view?
DoMatrixAllIter=True
# Separately visualize class averages of the last iteration?
DoShowLastIter=True
# Plot model (and mirror) fractions of the last iteration?
DoShowFractions=True
# Plot convergence statistics for all iterations?
DoShowStatsAllIter=True


# {hidden} Show expert options
"""If True, expert options will be displayed
"""
ShowExpertOptions = False

#------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
#------------------------------------------------------------------------------------------------

from protocol_dummy import *

if __name__ == '__main__':
    protocolMain(ProtDummy)
