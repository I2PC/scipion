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

#-----------------------------------------------------------------------------
# {section} Run parameters
#-----------------------------------------------------------------------------
# Run name:
""" This will identify your protocol run. It need to be unique for each protocol. You could have run1, run2 for protocol X, but not two
run1 for it. This name together with the protocol output folder will determine the working dir for this run.
"""
RunName = "run_001"

# {hidden} Run behavior
""" Resume from the last step, restart the whole process
"""
Behavior = "Resume"

#------------------------------------------------------------------------------------------
# {section} Picking parameters
#------------------------------------------------------------------------------------------
# {run}(preprocess_micrographs) Preprocessing Micrographs run
""" Directory with the preprocessing (output of the Preprocessing Micrographs protocol)
"""
PreprocessingRun = ""

# Perform automatic particle picking
""" Perform automatic particle picking """
AutomaticPicking = False

# {condition}(AutomaticPicking=True) Number of threads
""" This option provides shared-memory parallelization on multi-core machines.
It does not require any additional software, other than xmipp.
"""
NumberOfThreads = 2

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
