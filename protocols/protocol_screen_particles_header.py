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
# {file}(*.xmd *.stk *.sel *.ctfdat){validate}(PathExists) Input images
""" Stack file or metadata file with a selection file"""
InputFile=''

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