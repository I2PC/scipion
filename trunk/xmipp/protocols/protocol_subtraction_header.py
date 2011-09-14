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
#-----------------------------------------------------------------------------
# {section} Global parameters
#-----------------------------------------------------------------------------
#Comment
Comment='Describe your project here...'

#Comment
""" Subtraction subdirectory. 'run_001' will create 'subtraction/run_001' directory
"""
#RunName='run_008_exp_small'
RunName='run_009_exp_small'

# {list}(Resume, Restart) Run behavior
""" Resume from the last step, restart the whole process or continue at a given step or iteration
"""
Behavior = "Restart"
#Behavior = "Resume"
#Behavior = "Continue"

# {file} Protocol Name
ProtocolName='ProjMatch/xmipp2.4_exp_small/xmipp_protocol_projmatch_backup_test.py'
#ProtocolName='ProjMatch/xmipp_2.4_subtraction_crunchy/xmipp_protocol_projmatch_backup.py'

#Show results for iteration
""" Use data coming from iteration
"""
iterationNo=1

# {expert} Resume at Iter (vs Step)
"""This option control how to resume a previously performed run.
    Set to TRUE to restart at the beginning of iteration N
    or FALSE to continue at step N. (N is set in the next parameter).
    NOTE:If you do not know what are you doing make it equal to False
"""
IsIter =False

# Resume at iteration
""" Set to 1 to start a new run, set to -1 to continue the process (where you left it),
    set to a positive number N to restart at the begining of iteration N
    Note1: Do NOT delete working directory if this option is not set to 1
    Note2: Set this option to -1 if you want to perform extra iterations after
           successfully finish an execution
"""
ContinueAtIteration =1


#{expert}{file} Select the input volume
"""{expert DELETE THIS ENTRY}Name of the reference volume by default Iter_X_reconstruction.vol

"""
szInputVolumeName=''
#{expert}{file} Select docfile used to compute references
"""Name of the doc file used to compute reference library, by dfault
   ../Iter_(X-1)/Iter_(X-1)_current_angles.doc
"""
DocFileExp=''

# {expert} Mask reference volume?
doMask =True

#Crown Mask radius (inner)
dRradiusMin=39

#{expert}  iCrown mask radius center (outter)
dRradiusMax=64

# {expert} Backup Experimental Proj. angles
doBackupProjectionAngles =True

# {expert}Angular sampling rate
"""Angular distance (in degrees) between neighboring projection  points
   usually obtained from protocol file
"""   
AngSamplingRateDeg=''

# {expert}  Angular search range 
"""Maximum change in rot & tilt  (in +/- degrees)
   usually obtained from protocol file
"""   
MaxChangeInAngles =''

# {expert} Symmetry group
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry
    for a description of the symmetry groups format
    If no symmetry is present, give c1
    usually obtained from protocol file
"""
SymmetryGroup=''

CTFDatName=''
# {expert} Correct by CTF:
""" Set to True if you want to correct by CTF
"""
doCTFCorrection=False

# Scale images?
""" Set to True if you want to scale images (Using padding/windowing in Fourier space)
"""
doScaleImages=False

# New X dimension
""" New X dimension.
"""
dimX = 32

# {expert} New Y dimensions
""" New Y dimension. -1 means Y = X
"""
dimY = -1

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Number of (shared-memory) threads?
""" This option provides shared-memory parallelization on multi-core machines. 
    It does not require any additional software, other than xmipp
"""
NumberOfThreads = 2

# Number of MPI processes to use
NumberOfMpi = 2

#MPI job size 
"""Minimum size of jobs in mpi processes. Set to 1 for large images (e.g. 500x500) and to 10 for small images (e.g. 100x100)
"""
MpiJobSize ='1'

#Submit to queue
"""Submmit to queue"""
SubmitToQueue=False

#------------------------------------------------------------------------------------------------
# {section}{expert}{condition}(SubmitToQueue=True) Queue 
#------------------------------------------------------------------------------------------------

# Queue name
"""Name of the queue to submit the job"""
QueueName="high"

# Queue hours
"""This establish a maximum number of hours the job will
be running, after that time it will be killed by the 
queue system"""
QueueHours=72

#------------------------------------------------------------------------------------------------
# {hidden} Analysis of results

""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript='visualize_partial_projection_subtraction.py'

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
        
