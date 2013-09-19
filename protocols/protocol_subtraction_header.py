#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for subtraction
#
# Example use:
# ./protocol_subtraction_header.py
#
# Authors: Roberto Marabini,
#          Alejandro Echeverria Rey.
#
# {begin_of_header}

# {eval} expandCommentRun(allowContinue=True)

#-----------------------------------------------------------------------------
# {section} Input
#-----------------------------------------------------------------------------

# {run}(projmatch){post_run}(postSelectRunProjmatch) Import Projection Matching Run
"""
Select desired RUN from which you obtained the micrographs.

Possible input protocols are:
<Projection Matching>

"""
ImportRun = ''

#Use results for iteration
""" Use data coming from projection matching 
iteration number... Will be set to the last
 iteration each time a  projection matching
run is selected
"""
iterationNo=-1

#Use Reference Number
""" Use 3D reference from projection matching 
"""
reference3DNo=1

#{hidden}ReferenceFileName
"""{hiden}Name of 
the reference volume by default 
Iter_X_reconstruction.vol
"""
ReferenceFileNames=''

#{hidden} Select docfile used to compute references
"""Name of the doc file used to compute reference library, by default
   ../Iter_(X-1)/Iter_(X-1)_current_angles.doc
"""
DocFileExp=''

# {hidden} Mask reference volume?
"""mask reference volume? by default
use the same than in projection matching"""
doMask=True

# {expert}{condition}(doMask){wizard}(wizardSetAlignRadii) Crown mask radius center inner
"""mask reference volume? by default
use the same than in projection matching"""
InnerRadius=0

# {expert}{condition}(doMask){wizard}(wizardSetAlignRadii) Crown mask radius center outter
"""mask reference volume? by default
use the same than in projection matching"""
OuterRadius=-1

# {expert}Angular sampling rate
"""Angular distance (in degrees) between neighboring projection  points 
mask reference volume? by default
use the same than in projection matching
"""
AngSamplingRateDeg=-1

# {expert}  Angular search range 
"""Maximum change in rot & tilt  (in +/- degrees)
   usually obtained from projection matching protocol
"""   
MaxChangeInAngles =-1

# {expert} Symmetry group
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry
    for a description of the symmetry groups format
    If no symmetry is present, give c1
    usually obtained from protocol file
"""
SymmetryGroup=''

# {hidden} CTFDatInformation
"""mapping between particles and CTFs
"""
CTFDatName=''
# Correct by CTF:
""" Set to True if you want to correct by CTF
"""
doCTFCorrection=True

#-----------------------------------------------------------------------------
# {section}{has_question}  Scale images REMOVE?
#-----------------------------------------------------------------------------
#  Scale images
""" Scale image. Not supported yet..
"""
doScaleImages=True

# New X dimension
""" New X dimension.
"""
dimX = 64

# {expert} New Y dimensions
""" New Y dimension. -1 means Y = X
"""
dimY = -1

# {eval} expandParallel(mpi=3, jobsize=1)

#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript ='visualize_subtraction.py'

#------------------------------------------------------------------------------------------------
# {section}{visualize} Visualization
#------------------------------------------------------------------------------------------------

# {view} Display stack with reference images
""" Stack with reference images
"""
DisplayReference=False

# {view} Display stack with experimental images
""" Stack with experimental images
"""
DisplayExperimental=False

# {view} Display stack with subtracted images
""" Stack with subtracted images
"""
DisplaySubtracted=True

#-----------------------------------------------------------------------------
# {section} Debug
#-----------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

#Verify
"""Check that some output files are created. 
"""
Verify=True

# {expert} print wrapper name
PrintWrapperCommand=True

# {expert} print wrapper parameters
PrintWrapperParameters=True

# {expert} show file verification
ViewVerifyedFiles=True 

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
#SystemFlavour = "TORQUE-OPENMPI"
from protocol_subtraction import *
           
if __name__ == '__main__':
    protocolMain(ProtPartialProjectionSubtraction)
        
