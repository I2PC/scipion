#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based ML3D/MLF3D classification, according to:
#
# Example use:
# ./xmipp_protocol_ml3d.py
#
# Author: Sjors Scheres, December 2007
# Updated:  J. M. de la Rosa Trevin Sept 2011
#
#------------------------------------------------------------------------------------------------
# {begin_of_header}

#   {please_cite}
"""
for ML3D:  Scheres et al. (2007) Nature Methods, 4, 27-29
for MLF3D: Scheres et al. (2007) Structure, 15, 1167-1177
"""

#{include} inc_comment_run.py

#------------------------------------------------------------------------------------------
# {section} Input
#------------------------------------------------------------------------------------------
# {file}{validate}(PathExists) Input images:
""" 
Provide a list of images from a stack or metadata file that make up your data set.
The filenames should be relative to the <ProjectDir> where you are running the <Protocols>
If you have images outside the <Project> you should import them first.
"""
ImgMd = "all_images.xmd"

# {file}{validate}(PathExists) Initial 3D reference volumes:
""" Initial 3D density maps with the same dimensions as your particles.
"""
RefMd = "result_classes.xmd"

# Number of seeds per reference
'''
The total number of seeds generated will be the number of provided
reference volumes times the number of seeds per reference.
If you provide 2 initial volumes and 3 seeds per referece you will
produce 6 3D maps.
'''
NumberOfReferences = 1

#------------------------------------------------------------------------------------------------
# {section}{has_question} MLF parameters
#------------------------------------------------------------------------------------------------

# Perform MLF3D instead of ML3D
DoMlf = False

# Use CTF-amplitude correction inside MLF?
""" If set to true, provide the ctfdat file in the field below. 
If set to false, the ctfdat field will be ignored.
"""
DoCorrectAmplitudes = True

# {file}{validate}(PathExists){condition}(DoCorrectAmplitudes) CTFdat file with the input images:
""" The names of both the images and the ctf-parameter files should be 
relative path from <ProjectDir>. You can import data from other projects.
"""
InCtfDatFile='all_images.ctfdat'

# High-resolution limit (in Angstroms)
""" No frequencies higher than this limit will be taken into account. 
If zero is given, no limit is imposed.
"""
HighResLimit = 20

# Are the images CTF phase flipped?
""" You can run MLF with or without having phase flipped the images.
"""
ImagesArePhaseFlipped = True

# Are initial references CTF-amplitude corrected?
""" If coming from programs other than <xmipp_mlf_refine3d> this is 
usually not the case. If you will perform a grey-scale correction, 
this parameter becomes irrelevant as the output maps never have the
 CTF-amplitudes corrected.
"""
InitialMapIsAmplitudeCorrected = False

# {expert} Are the seeds CTF-amplitude corrected?
""" This option is only relevant if you provide your own seeds! 
If the seeds are generated automatically, this parameter becomes 
irrelevant as they will always be amplitude-corrected
"""
SeedsAreAmplitudeCorrected = False

# Correct the absolute grey-scale of initial references?
""" The probabilities are based on squared differences, so that the 
absolute grey scale is important.
"""
DoCorrectGreyScale = False
# {expert}{condition}(DoCorrectGreyScale) Samplig for projection matching
""" 
Angular sampling for a quick projection matching to obtain right grey scale.
As the resolution of the intial reference should be low, this sampling can
 be relatively crude, e.g. 15
"""
ProjMatchSampling = 15

# Low-pass filter initial references?
""" It is highly recommended to low-pass filter your initial reference 
volume as much as you can.
"""
DoLowPassFilterReference = True
# {condition}(DoLowPassFilterReference) Resolution of the low-pass filter (in Angstroms):
LowPassFilter = 50
# {condition}(DoLowPassFilterReference) Pixel size (in Angstroms):
PixelSize = 5.6

#------------------------------------------------------------------------------------------------
# {section}{has_question} ML3D classification
#------------------------------------------------------------------------------------------------
# Perform ML(F)3D classification
DoML3DClassification=True

# Angular sampling for classification:
""" 
Fine samplings take huge amounts of CPU and memory.
Therefore, in general, dont use samplings finer than 10 degrees.
"""
AngularSampling = 10

# Number of ML(F)3D iterations to perform:
NumberOfIterations = 25

# Point group symmetry:
""" 
See [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry]
for a description of the point group symmetries
Give <c1> if no symmetry is present
"""
Symmetry = 'c1'
# Refine the normalization for each image?
""" This variant of the algorithm deals with normalization errors.
    For more info see (and please cite):
    <Scheres et. al. (2009) J. Struc. Biol., in press>
"""
DoNorm = False

# {expert} Restart after iteration:
""" For previous runs that stopped before convergence,
    resume the calculations after the completely finished iteration,
    i.e. including all 3D reconstructions.
    Note that all flags about grey-scale correction, filtering and
    seed generation will be ignored if a value larger than 0 is given,
    since this option only concerns the ML3D classification part
"""
RestartIter=0
# {expert} Additional parameters:
""" 
Additional <xmipp_ml(f)_refine3d parameters>
For a complete description see the manual pages:
[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/MLrefine3D]
"""
ExtraParams=''

#-----------------------------------------------------------------------------
# {section}{expert} 3D Reconstruction
#-----------------------------------------------------------------------------

# {list}(wslART, fourier) Reconstruction method
""" Choose between wslART or fourier
"""
ReconstructionMethod ='wslART'

# {expert}{condition}(ReconstructionMethod=="wslART") Extra parameters
""" Additional reconstruction parameters for ART
    For details see:
    [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Art]
"""
ARTExtraParams =''

# {condition}(ReconstructionMethod=="fourier") Extra parameters
""" The Fourier-interpolation reconstruction method is much faster than wlsART 
    and may give similar results. It however is not guaranteed to optimize the 
    likelihood function. This is an experimental feature. One may limit the 
    maximum resolution of the fourier-interpolation using "-max_resolution 0.3"
    (to 0.3 digital frequency). Use the extra parameter entry below for that. 
    For details see:
    [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Fourier]
"""
FourierExtraParams =''

# {include} inc_parallel.py

# {hidden} Show expert options
"""If True, expert options will be displayed
"""
ShowExpertOptions = False

#------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
#------------------------------------------------------------------------------------------------

from protocol_ml3d import *

if __name__ == '__main__':
    protocolMain(ProtML3D)

