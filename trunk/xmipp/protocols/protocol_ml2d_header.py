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

#------------------------------------------------------------------------------------------
# {section} Input
#------------------------------------------------------------------------------------------
# {file}{validate}(PathExists) Input images metadata:
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR.
"""
ImgMd = "results_images.xmd"

# Generate references (or classes) ?
""" If you set to No, you should provide a metadata with the references images. 
The default generation is done by averaging subsets of input images
"""
DoGenerateReferences = True

# {condition}(DoGenerateReferences) Number of references:
""" Number of references to be generated. """
NumberOfReferences = 3

# {file}{validate}(PathExists){condition}(not DoGenerateReferences) References metadata:
""" Number of references to be generated. """
RefMd = "result_classes.xmd"

#------------------------------------------------------------------------------------------
# {section}{has_question} MLF-specific parameters
#------------------------------------------------------------------------------------------
# Use MLF2D instead of ML2D
DoMlf = False

# Use CTF-amplitude correction inside MLF?
""" If set to true, provide the ctfdat file in the field below. If set to false, one can ignore the ctfdat field, but has to provide the image pixel size in Angstrom
"""
DoCorrectAmplitudes = True

# {file}{condition}(DoCorrectAmplitudes) CTFdat file with the input images:
""" The names of both the images and the ctf-parameter files should be with absolute paths.
"""
InCtfDatFile = "all_images.ctfdat"

# {condition}(DoCorrectAmplitudes)Image pixel size (in Angstroms)
PixelSize = 4

# Are the images CTF phase flipped?
""" You can run MLF with or without having phase flipped the images.
"""
ImagesArePhaseFlipped = True

# High-resolution limit (in Angstroms)
""" No frequencies higher than this limit will be taken into account. If zero is given, no limit is imposed
"""
HighResLimit = 20

#------------------------------------------------------------------------------------------
# {section}{has_question} Advanced parameters
#------------------------------------------------------------------------------------------
# Show advanced parameters
AdvancedParameters = False

# Also include mirror in the alignment?
"""  Including the mirror transformation is useful if your particles have a handedness and may fall either face-up or face-down on the grid.
Note that when you want to use this ML2D run for later RCT reconstruction, you can NOT include the mirror transformation here.
"""
DoMirror = False

# {condition}(not DoMlf) Use the fast version of this algorithm?
""" See Scheres et al., Bioinformatics, 21 (Suppl. 2), ii243-ii244:
http://dx.doi.org/10.1093/bioinformatics/bti1140
"""
DoFast = True

# Refine the normalization for each image?
""" This variant of the algorithm deals with normalization errors. For more info see (and please cite) Scheres et. al. (2009) J. Struc. Biol., in press
"""
DoNorm = False

# Restart after iteration:
""" For previous runs that stopped before convergence, resume the calculations
after the completely finished iteration. (Use zero to start from the beginning)
"""
RestartIter = 0

# {expert} --eps   float
"""
 Stopping criterium 
"""
Eps = 5e-5

# {expert} --iter   int
"""
 Maximum number of iterations to perform  
"""
MaxIters = 100

#{expert}  --psi_step   float
"""
 In-plane rotation sampling interval [deg] 
"""
PsiStep = 5.0

#{expert}  --noise   float
"""
 Expected standard deviation for pixel noise  
"""
StdNoise = 1

# {expert} --offset   float
"""
 Expected standard deviation for origin offset [pix] 
"""
StdOffset = 3.0

# {expert} -C   double
"""
 Significance criterion for fast approach  
"""
CFast = 1e-12

# {expert} --zero_offsets
"""
 Kick-start the fast algorithm from all-zero offsets  
"""
ZeroOffsets = False

# {expert} --fix_sigma_noise
"""
 Do not re-estimate the standard deviation in the pixel noise  
"""
FixSigmaNoise = False
# {expert}--fix_sigma_offset
"""
 Do not re-estimate the standard deviation in the origin offsets  
"""
FixSigmaOffset = False
# {expert} --fix_fractions
"""
 Do not re-estimate the model fractions  
"""
FixFractions = False

#------------------------------------------------------------------------------------------
# {section} Parallelization 
#------------------------------------------------------------------------------------------
#  {condition}(not DoMlf) Number of threads
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

# {condition}(SubmitToQueue) Queue hours
"""This establish a maximum number of hours the job will
be running, after that time it will be killed by the
queue system
"""
QueueHours = 72

#------------------------------------------------------------------------------------------------
# {section}{visualize} Visualization
#------------------------------------------------------------------------------------------------
# Visualize the class averages of all iterations in matrix-view?
#DoMatrixAllIter = True
# Separately visualize class averages of the last iteration?
#DoShowLastIter = True
# Plot model (and mirror) fractions of the last iteration?
#DoShowFractions = True
# Plot convergence statistics for all iterations?
#DoShowStatsAllIter = True


# {hidden} Show expert options
"""If True, expert options will be displayed
"""
ShowExpertOptions = False

#------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
#------------------------------------------------------------------------------------------------

from protocol_ml2d import *

if __name__ == '__main__':
    protocolMain(ProtML2D)
