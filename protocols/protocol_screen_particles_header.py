#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based screening single-particles: 
# Author: Carlos Oscar, October 2011
#
# {begin_of_header}

#{include} inc_comment_run.py

# {eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Screening parameters
#-----------------------------------------------------------------------------
# {file}(images*.xmd){validate}(PathExists) Input images
""" Stack file or metadata file with a selection file"""
InputFile=''

# {expert}{list_combo}(None, MaxZscore, Percentage) Automatic particle rejection
""" How to automatically reject particles. It can be none (no rejection), maxZscore (reject a particle
if its Zscore is larger than this value), Percentage (reject a given percentage in each one of the screening criteria). """
RejectionMethod='None'

# {expert}{condition}(RejectionMethod=="MaxZscore") Maximum Zscore
MaxZscore=3

# {expert}{condition}(RejectionMethod=="Percentage") Percentage (%)
Percentage=5

# {eval} expandParallel(mpi=0, threads=0, hours=4)
# {eval} expandExpert()

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
