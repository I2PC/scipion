#!/usr/bin/env xmipp_python
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

#------------------------------------------------------------------------------------------
# {begin_of_header}

#{eval} expandCommentRun()

#  {cite} ML References
CiteML = """
for ML2D:  Scheres et al. (2005) J.Mol.Biol 348, 139-149
for MLF2D: Scheres et al. (2007) Structure 15, 1167-1177
"""

#------------------------------------------------------------------------------------------
# {section} Input
#------------------------------------------------------------------------------------------

# {file}(images*.xmd){validate}(PathExists) Input images:
""" 
Provide a list of images from a stack or metadata file that make up 
your data set. The filenames should be relative to the <ProjectDir> 
where you are running the <Protocols>
If you have images outside the <Project> you should import them first.
"""
ImgMd = ""

# Generate references (or classes) ?
""" 
If you set to <No>, you should provide references images 
(In stack or metadata file). 
The default generation is done by averaging subsets 
of input images.
"""
DoGenerateReferences = True

# {condition}(DoGenerateReferences) Number of references:
""" Number of references to be generated. """
NumberOfReferences = 3

# {file}(result_classes*xmd){validate}(PathExists){condition}(not DoGenerateReferences) References image(s):
""" Metadata with the input reference images """
RefMd = ""

#------------------------------------------------------------------------------------------
# {section}{has_question} MLF-specific parameters
#------------------------------------------------------------------------------------------
# Use MLF2D instead of ML2D
DoMlf = False

# Use CTF-amplitude correction inside MLF?
""" 
If set to <Yes>, the input images file should contains
the CTF information for each image.
If set to <No>, provide the images pixel size in Angstrom
"""
DoCorrectAmplitudes = True

# Are the images CTF phase flipped?
""" You can run MLF with or without having phase flipped the images.
"""
ImagesArePhaseFlipped = True

# High-resolution limit (in Angstroms)
""" 
No frequencies higher than this limit will be taken into account. 
If zero is given, no limit is imposed
"""
HighResLimit = 20

#------------------------------------------------------------------------------------------
# {section}{has_question} Advanced parameters
#------------------------------------------------------------------------------------------
# Show advanced parameters
AdvancedParameters = False

# Also include mirror in the alignment?
"""
Including the mirror transformation is useful if your particles 
have a handedness and may fall either face-up or face-down on the grid. 
Note that when you want to use this <ML2D> run for later <RCT> 
reconstruction, you can <NOT include the mirror> transformation here.
"""
DoMirror = True

# {condition}(not DoMlf) Use the fast version of this algorithm?
"""
For details see:
<See Scheres et al., Bioinformatics, 21 (Suppl. 2), ii243-ii244>
[http://dx.doi.org/10.1093/bioinformatics/bti1140]
"""
DoFast = True

# Refine the normalization for each image?
""" 
This variant of the algorithm deals with normalization errors. 
For more info see (and please cite):
<Scheres et. al. (2009) J. Struc. Biol., in press>
"""
DoNorm = False

# Restart after iteration:
""" 
For previous runs that stopped before convergence, resume the 
calculations after the completely finished iteration. 
(Use <zero> to start from the beginning)
"""
RestartIter = 0

# {expert} Stoping threshold
"""
When the differences beetween the model on the previous iteration
and the current one is less than this threshold, it will be said
that convergence has been reached
"""
Eps = 5e-5

# {expert} Maximum number of iterations
"""
If the convergence has not been reached after this number 
of iterations, the process will be stopped.
"""
MaxIters = 100

#{expert}(PsiStep >= 1.0)  In-plane rotation sampling (degrees)
"""
In-plane rotation sampling interval (degrees) 
"""
PsiStep = 5.0

#{expert}  Std for pixel noise
"""
Expected standard deviation for pixel noise  
"""
StdNoise = 1.0

# {expert} Std for origin offset
"""
Expected standard deviation for origin offset (pixels) 
"""
StdOffset = 3.0


#{eval} expandParallel()

#------------------------------------------------------------------------------------------------
# {section}{visualize} Visualization
#------------------------------------------------------------------------------------------------
# {view} Visualize last iter references
DoShowReferences = True

# {view} Show Log-Likelihood over iterations?
""" The Log-Likelihood value should increase"""
DoShowLL = True

# {view} Show maximum model probability?
""" 
Show the maximum probability for a model, this 
should tend to be a deltha function
"""
DoShowPmax = True

# {view} Show plot for signal change?
"""Should approach to zero when convergence"""
DoShowSignalChange = True

# {view} {condition}(DoMirror) Show mirror fraction of last iteration?
"""The the mirror fraction of each reference in last iteration"""
DoShowMirror = True

#------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
#------------------------------------------------------------------------------------------------

from protocol_ml2d import *

if __name__ == '__main__':
    protocolMain(ProtML2D)
