#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# {begin_of_header}

#{eval} expandCommentRun()

#------------------------------------------------------------------------------------------------
# {section} Input  
#------------------------------------------------------------------------------------------------
# {file}(*.emx) EMX file name
"""EMX file name:"""
EmxFileName = ''

# {validate}(IsFloatOrEmpty) Sampling rate (A/pixel)
""" If left blank, the value will be read from the EMX file. """
SamplingRate = ""

# Sort particles?
"""
If set to Yes, the particles will be sorted by similarity.
A ZScore column will be added to output <images.xmd>.
A higher ZScore means a more different particles.
""" 
DoSort = False

# {eval} expandParticlesPreprocess(allowFlip=False)

# {eval} expandParallel(threads=0, mpi=8, hours=12)

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_emx_import_particles import *

if __name__ == '__main__':
    protocolMain(ProtEmxImportParticles)
