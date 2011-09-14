#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using hierarchical clustering principles
#
# Example use:
# ./xmipp_protocol_cl2d.py
#
# Author: Carlos Oscar Sanchez Sorzano, July 2009
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
# {section} CL2D parameters
#------------------------------------------------------------------------------------------------
# {file}{validate}(PathExists) Selfile with the input images:
""" This selfile points to the stack or metadata containing your images 
"""
InSelFile=''

# Number of references (or classes) to be used:
NumberOfReferences=64

# {expert} Number of initial references
""" Initial number of initial models
"""
NumberOfInitialReferences=4

# {expert} Number of iterations
""" Maximum number of iterations within each level
"""
NumberOfIterations=15

# {expert}{list_combo}(correlation, correntropy) Comparison method
""" Use correlation or correntropy """
ComparisonMethod='correlation'

# {expert}{list_combo}(classical, robust) Clustering criterion
""" Use the classical clustering criterion or the robust clustering criterion """
ClusteringMethod='classical'

# {expert} Additional parameters for classify_CL2D
""" --verbose, --corrSplit, ...
"""
AdditionalParameters=''

#------------------------------------------------------------------------------------------------
# {section}{expert}{condition}(NumberOfReferences>NumberOfInitialReferences) Core analysis
#------------------------------------------------------------------------------------------------
# Junk Zscore
""" Which is the average Z-score to be considered as junk. Typical values
    go from 1.5 to 3. For the Gaussian distribution 99.5% of the data is
    within a Z-score of 3. Lower Z-scores reject more images. Higher Z-scores
    accept more images."""
thZscore=3

# PCA Zscore
""" Which is the PCA Z-score to be considered as junk. Typical values
    go from 1.5 to 3. For the Gaussian distribution 99.5% of the data is
    within a Z-score of 3. Lower Z-scores reject more images. Higher Z-scores
    accept more images.
    
    This Z-score is measured after projecting onto the PCA space."""
thPCAZscore=3

# Tolerance
""" An image belongs to the stable core if it has been with other images in the same class
in all the previous levels except possibly a few of them. Tolerance defines how few is few.
Tolerance=0 means that an image must be in all previous levels with the rest of images in
the core."""
Tolerance=1

#------------------------------------------------------------------------------------------
# {section} Parallelization
#------------------------------------------------------------------------------------------
# Number of MPI processes
""" Set to 1 if you do not have MPI installed"""
NumberOfMpi = 4

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

#------------------------------------------------------------------------------------------------
# {section}{visualize} Visualization
#------------------------------------------------------------------------------------------------
# {list_combo}(Classes,Class Cores,Class Stable Cores) What to show
WhatToShow="Class Cores"

# Visualize last level
DoShowLast=True

# {condition}(DoShowLast==False) Visualize several levels
"""Create a  list of levels like: 0,1,3 or 0-3 """
LevelsToShow=""

# {hidden} Show expert options
"""If True, expert options will be displayed
"""
ShowExpertOptions = False

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#		
# Main
#     
from protocol_cl2d import *

if __name__ == '__main__':
    protocolMain(ProtCL2D)
