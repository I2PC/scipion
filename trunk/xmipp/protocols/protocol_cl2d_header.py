#!/usr/bin/env python
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

# {file} Selfile with the input images:
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
InSelFile='sort_junk.sel'


#------------------------------------------------------------------------------------------------
# {section} Class averages parameters
#------------------------------------------------------------------------------------------------
# Number of references (or classes) to be used:
NumberOfReferences=64

# {expert} Number of initial references
""" Initial number of initial models
"""
NumberOfReferences0=4

# {expert} Number of iterations
""" Maximum number of iterations within each level
"""
NumberOfIterations=15

# {expert} Band pass filter
""" Apply a band pass filter before clustering """
DoFilter =True

# {expert}{condition}(DoFilter=True) Highpass cutoff frequency
""" In (Angstroms/Pixel). Set to 0 if not desired """
Highpass =0.02

# {expert}{condition}(DoFilter=True) Lowpass cutoff frequency
""" In (Angstroms/Pixel). Set to 0 if not desired """
Lowpass =0.4

# {expert}{list}(correlation, correntropy) Comparison method
""" Use correlation or correntropy """
ComparisonMethod='correlation'

# {expert}{list}(classical, robust) Clustering criterion
""" Use the classical clustering criterion or the robust clustering criterion """
ClusteringMethod='classical'

# {expert} Additional parameters for classify_CL2D
""" -verbose, -corrSplit, ...
"""
AdditionalParameters=''

#------------------------------------------------------------------------------------------------
# {section} Core analysis
#------------------------------------------------------------------------------------------------
# Good class core size (%)
""" A class is a good class if at least this percentage (around 50%) of the
    images assigned to it have been together in all the previous levels.
    Larger values of this parameter tend to keep few good classes. Smaller
    values of this parameter tend to consider more classes as good ones."""
thGoodClass=50

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

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#		
# Main
#     
from protocol_cl2d import *

if __name__ == '__main__':
    protocolMain(ProtCL2D)
