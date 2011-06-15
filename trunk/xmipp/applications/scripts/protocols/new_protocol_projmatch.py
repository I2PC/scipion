#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for projection matching
#
# Example use:
# ./xmipp_protocol_projmatch.py
#
# Authors: Roberto Marabini,
#          Sjors Scheres,    March 2008
#        Rewritten by Roberto Marabini
#-----------------------------------------------------------------------------
# {section} Global parameters
#-----------------------------------------------------------------------------
#Comment
Comment='Describe your project here...'
# {file} Selfile with the input images:
#from XmippData import SingleImgSize
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
SelFileName ='new20.sel'

# {file} {expert} Docfile with the input angles:
""" Do not provide anything if there are no angles yet. 
    In that case, the starting angles will be read from the image headers
    This docfile should be in newXmipp-style format (with filenames as comments)
    Note that all filenames in this docfile should be with absolute paths!
"""
DocFileName =' '

# {file} Initial 3D reference map:
""" Write down the reference/es name. For example "Reference1.vol Reference2.vol"
    specifies two references
"""
ReferenceFileNames ='ico1.vol ico2.vol ico3.vol'

# Working subdirectory: 
""" This directory will be created if it doesn't exist, and will be used to store all output from this run. Don't use the same directory for multiple different runs, instead use a structure like run1, run2 etc. 
"""
WorkingDir ='ProjMatch/new20'

# Delete working subdirectory if it already exists?
""" Just be careful with this option...
"""
DoDeleteWorkingDir =True

# Number of iterations to perform
NumberofIterations = 4

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

# {expert} Save disc space by cleaning up intermediate files?
""" Be careful, many options of the visualization protocol will not work anymore, 
    since all class averages, selfiles etc will be deleted.
"""
CleanUpFiles =False

# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project. Often, each data set of a given sample has its own ProjectDir.
ProjectDir='/gpfs/fs1/home/bioinfo/roberto/PhantomIco'
"""
ProjectDir='/gpfs/fs1/home/bioinfo/roberto/PhantomIco'

# {expert} Directory name for logfiles:
LogDir ='Logs'

#-----------------------------------------------------------------------------
# {section} CTF correction
#-----------------------------------------------------------------------------

# Perform CTF correction?
""" If set to true, a CTF (amplitude and phase) corrected map will be refined,
    and the data will be processed in CTF groups.
    Note that you cannot combine CTF-correction with re-alignment of the classes.
"""
DoCtfCorrection =True

# {file} CTFDat file with CTF data:
""" The input selfile may be a subset of the images in the CTFDat file, but all 
    images in the input selfile must be present in the CTFDat file. This field is 
    obligatory if CTF correction is to be performed. 
    Note that this file should be positioned in the project directory, and that the
    image names and ctf parameter filenames should be in absolute paths.
"""
CTFDatName ='new_ctf.ctfdat'

# Make CTF groups automatically?
""" Make CTF groups based on a maximum differences at a given resolution limit.
    If this option is set to false, a docfile with the defocus values where to 
    split the images in distinct defocus group has to be provided (see expert option below)
"""
DoAutoCtfGroup =True

# Maximum difference in CTF-values in one group
""" If the difference between the CTF-values up to the resolution limit specified 
    below is larger than the value given here, two images will be placed in 
    distinct CTF groups.
"""
CtfGroupMaxDiff = 0.1

# Resolution limit (Ang) for CTF-differences in one group
""" Maximum resolution where to consider CTF-differences among different groups.
    One should use somewhat higher resolutions than those aimed for in the refinement.
"""
CtfGroupMaxResol = 5.6

# {file} {expert} Docfile with defocus values where to split into groups
""" This field is obligatory if you do not want to make the CTF groups automatically.
    Note that the requested docfile can be made initially with the xmipp_ctf_group program,
    and then it can be edited manually to suit your needs. 
"""
SplitDefocusDocFile =''

# {expert} Padding factor
""" Application of CTFs to reference projections and of Wiener filter to class averages will be done using padded images.
    Use values larger than one to pad the images. Suggestion, use 1 for large image and 2 for small
"""
PaddingFactor =2

# {expert} Wiener constant
""" Term that will be added to the denominator of the Wiener filter.
    In theory, this value is the inverse of the signal-to-noise ratio
    If a negative value is taken, the program will use a default value as in FREALIGN 
    (i.e. 10% of average sum terms over entire space) 
    see Grigorieff JSB 157 (2006) pp117-125
"""
WienerConstant = -1

# Images have been phase flipped?
DataArePhaseFlipped =True

# Is the initial reference map CTF (amplitude) corrected?
"""
    You may specify this option for each iteration. 
    This can be done by a sequence of 0 or 1 numbers (for instance, "1 1 0 0" 
    specifies 4 iterations, the first two applied alig2d while the last 2
    dont. an alternative compact notation is 
    is ("2x1 2x0", i.e.,
    2 iterations with value 1, and 2 with value 0).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
    IMPORTANT: if you set this variable to 0 the output  of the projection
    muching step will be copied as output of align2d
"""
ReferenceIsCtfCorrected ='1'

#-----------------------------------------------------------------------------
# {section} Mask
#-----------------------------------------------------------------------------
# Mask reference volume?
""" Masking the reference volume will increase the signal to noise ratio.
    Do not provide a very tight mask.
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mask for details
"""
DoMask =True

# Use a spherical mask?
""" If set to true, provide the radius of the mask in the next input field
    if set to false, provide a binary mask file in the second next input field
"""
DoSphericalMask =True

# Radius of spherical mask
""" This is the radius (in pixels) of the spherical mask 
"""
MaskRadius = 64

# {file} Binary mask file
""" This should be a binary (only 0/1-valued) Xmipp volume of equal dimension as your reference
    The protein region should be white (1) and the solvent should be black (0).
    Note that this entry is only relevant if no spherical mask is used.
"""
MaskFileName ='mask.vol'

#-----------------------------------------------------------------------------
# {section} Projection Matching
#-----------------------------------------------------------------------------
# Inner radius for rotational correlation:
""" In pixels from the image center
    You may specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "8 8 2 2 " 
    specifies 4 iterations, the first two set the value to 8 
    and the last two to 2. An alternative compact notation 
    is ("2x8 2x0", i.e.,
    2 iterations with value 8, and 2 with value 2).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
InnerRadius = '0'

# Outer radius for rotational correlation
""" In pixels from the image center. Use a negative number to use the entire image.
    WARNING: this radius will be use for masking before computing resolution
    You may specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "8 8 2 2 " 
    specifies 4 iterations, the first two set the value to 8 
    and the last two to 2. An alternative compact notation 
    is ("2x8 2x0", i.e.,
    2 iterations with value 8, and 2 with value 2).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
OuterRadius = '64'

# {expert} Available memory to store all references (Gb)
""" This is only for the storage of the references. If your projections do not fit in memory, 
    the projection matching program will run MUCH slower. But, keep in mind that probably 
    some additional memory is needed for the operating system etc.
    Note that the memory per computing node needs to be given. That is, when using threads, 
    this value will be multiplied automatically by the number of (shared-memory) threads.
"""
AvailableMemory = 2

# Angular sampling rate
"""Angular distance (in degrees) between neighboring projection  points
    You may specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "8 8 2 2 " 
    specifies 4 iterations, the first two set the value to 8 
    and the last two to 2. An alternative compact notation 
    is ("2x8 2x0", i.e.,
    2 iterations with value 8, and 2 with value 2).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
AngSamplingRateDeg='1 3 2 1'

# Angular search range 
"""Maximum change in rot & tilt  (in +/- degrees)
    You may specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "1000 1000 10 10 " 
    specifies 4 iterations, the first two set the value to 1000 (no restriction)
    and the last two to 10degrees. An alternative compact notation 
    is ("2x1000 2x10", i.e.,
    2 iterations with value 1000, and 2 with value 10).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
MaxChangeInAngles='1000 10 4 2'

# {expert} Perturb projection directions?
""" If set to 1, this option will result to a Gaussian perturbation to the 
    evenly sampled projection directions of the reference library. 
    This may serve to decrease the effects of model bias.
    You may specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "1 1 0" 
    specifies 3 iterations, the first two set the value to 1 
    and the last to 0. An alternative compact notation 
    is ("2x1 0", i.e.,
    2 iterations with value 1, and 1 with value 0).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
PerturbProjectionDirections ='0'

# Maximum change in origin offset
""" Maximum allowed change in shift in the 3D+2D searches (in +/- pixels).
    Shifts larger than this value will be reset to (0,0)
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "1000 1000 10 10 " 
    specifies 4 iterations, the first two set the value to 1000 (no restriction)
    and the last two to 10degrees. An alternative compact notation 
    is ("2x1000 2x10", i.e.,
    2 iterations with value 1000, and 2 with value 10).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
MaxChangeOffset='1000 10 5'

# Search range for 5D translational search 
""" Give search range from the image center for 5D searches (in +/- pixels).
    Values larger than 0 will results in 5D searches (which may be CPU-intensive)
    Give 0 for conventional 3D+2D searches. 
    Note that after the 5D search, for the optimal angles always 
    a 2D exhaustive search is performed anyway (making it ~5D+2D)
    Provide a sequence of numbers (for instance, "5 5 3 0" specifies 4 iterations,
    the first two set the value to 5, then one with 3, resp 0 pixels.
    An alternative compact notation is ("3x5 2x3 0", i.e.,
    3 iterations with value 5, and 2 with value 3 and the rest with 0).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
    
"""
Search5DShift ='4x5 0'

# {expert} Step size for 5D translational search
""" Provide a sequence of numbers (for instance, "2 2 1 1" specifies 4 iterations,
    the first two set the value to 2, then two with 1 pixel.
    An alternative compact notation is ("2x2 2x1", i.e.,
    2 iterations with value 2, and 2 with value 1).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
    
"""
Search5DStep ='2'

# {expert} Restrict tilt angle search?
DoRestricSearchbyTiltAngle =False

# {expert} Lower-value for restricted tilt angle search
Tilt0 = 40

# {expert} Higher-value for restricted tilt angle search
TiltF = 90

# Symmetry group
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry
    for a description of the symmetry groups format
    If no symmetry is present, give c1
"""
SymmetryGroup ='i3'

# {expert} Symmetry group for Neighbourhood computations
""" If you do not know what this is leave it blank.
    This symmetry will be using for compute neighboring points,
    but not for sampling or reconstruction
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry
    for a description of the symmetry groups format
    If no symmetry is present, give c1
"""
SymmetryGroupNeighbourhood =''

# {expert} compute only closest neighbor 
""" This option is only relevant if SymmetryGroupNeighbourhood !=''
    If set to 1 only one neighbor will be computed per sampling point
    You may specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "1 1 0" 
    specifies 3 iterations, the first two set the value to 1 
    and the last to 0. An alternative compact notation 
    is ("2x1 0", i.e.,
    2 iterations with value 1, and 1 with value 0).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
OnlyWinner ='0'

# Discard images with ccf below
""" Provide a sequence of numbers (for instance, "0.3 0.3 0.5 0.5" specifies 4 iterations,
    the first two set the value to 0.3, then two with 0.5.
    An alternative compact notation would be ("2x0.3 2x0.5").
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
    Set to -1 to prevent discarding any images
"""
MinimumCrossCorrelation ='-1'

# Discard percentage of images with ccf below
""" Provide a sequence of numbers (for instance, "20 20 10 10" specifies 4 iterations,
    the first two set the value to 20%, then two with 10%
    An alternative compact notation would be ("2x20 2x10").
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
    Set to zero to prevent discarding any images
"""
DiscardPercentage ='10'

# Perform scale search?
""" If true perform scale refinement
"""
DoScale =False

# Step scale factors size
""" Scale step factor size (1 means 0.01 in/de-crements arround 1)
"""
ScaleStep ='1'

# Number of scale steps
""" Number of scale steps.
    With default values (ScaleStep='1' and ScaleNumberOfSteps='3'): 1 +/-0.01 | +/-0.02 | +/-0.03.    
    With values ScaleStep='2' and ScaleNumberOfSteps='4' it performs a scale search over:
     1 +/-0.02 | +/-0.04 | +/-0.06 | +/-0.08.    
"""
ScaleNumberOfSteps ='3'


# {expert} Additional options for Projection_Matching
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Projection_matching and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_projection_matching for details
    try -Ri xx -Ro yy for restricting angular search (xx and yy are
    the particle inner and outter radius)
    
"""
ProjMatchingExtra =''

#-----------------------------------------------------------------------------
# {section} 2D re-alignment of classes
#-----------------------------------------------------------------------------
# Perform 2D re-alignment of classes?
""" After performing a 3D projection matching iteration, each of the
    subsets of images assigned to one of the library projections is
    re-aligned using a 2D-alignment protocol.
    This may serve to remove model bias.
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Align2d for details
    Note that you cannot combine this option with CTF-correction!
    You may specify this option for each iteration. 
    This can be done by a sequence of 0 or 1 numbers (for instance, "1 1 0 0" 
    specifies 4 iterations, the first two applied alig2d while the last 2
    dont. an alternative compact notation is 
    is ("2x1 2x0", i.e.,
    2 iterations with value 1, and 2 with value 0).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
    IMPORTANT: if you set this variable to 0 the output  of the projection
    muching step will be copied as output of align2d
"""
DoAlign2D ='0'

# {expert} Number of align2d iterations:
""" Use at least 3 iterations
    The number of align iteration may change in each projection matching iteration
    Ffor instance, "4 4 3 3 " 
    specifies 4 alig2d iterations in the first projection matching iteration 
    and  two 3 alig2d iteration in the last 2 projection matching iterations.
     An alternative compact notation 
    is ("2x4 2x3", i.e.,
    2 iterations with value 4, and 2 with value 3).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
Align2DIterNr ='4'

# {expert} Maximum change in origin offset (+/- pixels)
"""Maximum change in shift  (+/- pixels)
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "1000 1000 10 10 " 
    specifies 4 iterations, the first two set the value to 1000 (no restriction)
    and the last two to 10degrees. An alternative compact notation 
    is ("2x1000 2x10", i.e.,
    2 iterations with value 1000, and 2 with value 10).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
Align2dMaxChangeOffset ='2x1000 2x10'

# {expert} Maximum change in rotation (+/- degrees)
"""Maximum change in shift  (+/- pixels)
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "1000 1000 10 10 " 
    specifies 4 iterations, the first two set the value to 1000 (no restriction)
    and the last two to 10degrees. An alternative compact notation 
    is ("2x1000 2x10", i.e.,
    2 iterations with value 1000, and 2 with value 10).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
Align2dMaxChangeRot ='2x1000 2x20'

#-----------------------------------------------------------------------------
# {section} 3D Reconstruction
#-----------------------------------------------------------------------------

# {list}|fourier|art|wbp| Reconstruction method
""" Choose between wbp, art or fourier
"""
ReconstructionMethod ='fourier'

# {expert} Values of lambda for art
""" IMPORTANT: ou must specify a value of lambda for each iteration even
    if art has not been selected.
    IMPORTANT: NOte that we are using the WLS version of ART that 
    uses geater lambdas than the plain art.
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Art
        for details
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, ".1 .1 .3 .3" 
    specifies 4 iterations, the first two set the value to 0.1 
    (no restriction)
    and the last  two to .3. An alternative compact notation 
    is ("2x.1 2x.3").
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
ARTLambda ='0.2'

# {expert} Additional reconstruction parameters for ART
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Art
        for details
"""
ARTReconstructionExtraCommand ='-k 0.5 -n 10 '

# Initial maximum frequency used by reconstruct fourier
""" This number os only used in the first iteration. 
    From then on, it will be set to resolution computed in the resolution section
"""
FourierMaxFrequencyOfInterest =0.25

# {expert} Additional reconstruction parameters for WBP
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Wbp and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_wbp and
        for details
"""
WBPReconstructionExtraCommand =''

# {expert} Additional reconstruction parameters for Fourier
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Fourier and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_Fourier and
        for details
    -thr_width 
"""
FourierReconstructionExtraCommand =''

#-----------------------------------------------------------------------------
# {section} Compute Resolution
#-----------------------------------------------------------------------------
# Compute resolution?
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Resolution for details
    You may specify this option for each iteration. 
    This can be done by a sequence of 0 or 1 numbers (for instance, "1 1 0 0" 
    specifies 4 iterations, the first two applied alig2d while the last 2
    dont. an alternative compact notation is 
    is ("2x1 2x0", i.e.,
    2 iterations with value 1, and 2 with value 0).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
DoComputeResolution ='1'

# {expert} Split references averages
"""In theory each reference average should be splited
   in two when computing the resolution. In this way each
   projection direction will be represented in each of the
   subvolumes used to compute the resolution. A much faster
   but less accurate approach is to split the 
   proyection directions in two but not the averages. We
   recomend the first approach for small volumes and the second for
   large volumes (especially when using small angular
   sampling rates.
   IMPORTANT: the second option has ONLY been implemented for FOURIER
   reconstruction method. Other reconstruction methods require this
   flag to be set to True
    You may specify this option for each iteration. 
    This can be done by a sequence of 0 or 1 numbers (for instance, "1 1 0 0" 
    specifies 4 iterations, the first two split the images   while the last 2
    don't. an alternative compact notation is 
    is ("2x1 2x0", i.e.,
    2 iterations with value 1, and 2 with value 0).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more vapplications/scripts/protocols/new_protocol_projmatch.pyalues than iterations the extra value are ignored
"""
DoSplitReferenceImages ="1"


# Pixel size (in Ang.)
""" This will make that the X-axis in the resolution plots has units 1/Angstrom
"""
ResolSam=5.6

#-----------------------------------------------------------------------------
# {section} Postprocessing
#-----------------------------------------------------------------------------
# Low-pass filter the reference?
DoLowPassFilter =True

# Use estimated resolution for low-pass filtering?
"""If set to true, the volume will be filtered at a frecuency equal to
   the  resolution computed with a FSC=0.5 threshold, possibly 
   plus a constant provided by the user in the next input box. 

   If set to false, then the filtration will be made at the constant 
   value provided by the user in the next box (in digital frequency, 
   i.e. pixel-1: minimum 0, maximum 0.5) 
"""
UseFscForFilter =True

# Constant to by add to the estimated resolution
""" The meaning of this field depends on the previous flag.
    If set to true, then the volume will be filtered at a frecuency equal to
    the  resolution computed with resolution_fsc (FSC=0.5) plus the value 
    provided in this field 
    If set to false, the volume will be filtered at the resolution
    provided in this field 
    This value is in digital frequency, or pixel^-1: minimum 0, maximum 0.5

    You can specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, ".15 .15 .1 .1" 
    specifies 4 iterations, the first two set the constant to .15
    and the last two to 0.1. An alternative compact notation 
    is ("2x.15 2x0.1", i.e.,
    4 iterations with value 0.15, and three with value .1).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
ConstantToAddToFiltration ='0.1'

# {expert} Center volume
""" Center volume after each iteration """
DoCenterVolume =False

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Number of (shared-memory) threads?
""" This option provides shared-memory parallelization on multi-core machines. 
    It does not require any additional software, other than xmipp
"""
NumberOfThreads = 1

# distributed-memory parallelization (MPI)?
""" This option provides distributed-memory parallelization on multi-node machines. 
    It requires the installation of some MPI flavour, possibly together with a queueing system
"""
DoParallel =True

# Number of MPI processes to use:
NumberOfMpiProcesses =10

# minumum size of jobs in mpi processe. Set to 1 for large images (e.g. 500x500) and to 10 for small images (e.g. 100x100)
MpiJobSize ='1'

# MPI system Flavour 
""" Depending on your queuing system and your mpi implementation, different mpirun-like commands have to be given.
    Ask the person who installed your xmipp version, which option to use. 
    Or read: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/ParallelPage. The following values are available: 
"""
SystemFlavour ='TORQUE-OPENMPI'

#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript ='visualize_projmatch.py'
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

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#-----------------------------------------------------------------------------

from protocol_base import *

class ProtProjMatch(XmippProtocol):
    
    #Some class variables
    ReferenceVolumeName = 'reference_volume.vol'
    LibraryDir = "ReferenceLibrary"
    ProjectLibraryRootName = LibraryDir + "/gallery"
    ProjMatchDir = "ProjMatchClasses"
    ProjMatchName = 'proj_match'
    ClassAverageName = 'class_average'
    #ProjMatchRootName = ProjMatchDir + "/" + ProjMatchName
    ForReconstructionSel = "reconstruction.sel"
    ForReconstructionDoc = "reconstruction.doc"
    MultiAlign2dSel = "multi_align2d.sel"
    DocFileWithOriginalAngles = 'original_angles.doc'
    docfile_with_current_angles = 'current_angles'
    FilteredReconstruction = "filtered_reconstruction"
    
    ReconstructedVolume = "reconstruction"#
    maskReferenceVolume = "masked_reference"#
    
    OutputFsc = "resolution.fsc"
    CtfGroupDirectory = "CtfGroups"
    CtfGroupRootName = "ctf"
    CtfGroupSubsetFileName = CtfGroupRootName + "_images.sel"
    
    reconstructedFileNamesIters = []# names for reconstructed volumes
    #maskedFileNamesIter = []# names masked volumes used as reference
    numberOfReferences = 1#number of references
    createAuxTable = False
    NumberOfCtfGroups = 1
    
    def __init__(self, scriptname, workingdir, projectdir=None, restartStep=1, logdir='Logs'):
        super(Prot2,self).__init__(scriptname, workingdir, projectdir, restartStep, logdir)
    
    def preRun(self):
        #Convert directories/files  to absolute path from projdir
        self.CtfGroupDirectory = os.path.join(self.WorkingDir, self.CtfGroupDirectory)
        self.CtfGroupSubsetFileName = os.path.join(self.CtfGroupDirectory, self.CtfGroupSubsetFileName)
    #vector for iterations??????
    #    global ProjMatchDir
    #    ProjMatchDir = WorkingDir +'/' + ProjMatchDir
        self.DocFileWithOriginalAngles = os.path.join(self.WorkingDir, self.DocFileWithOriginalAngles)
    
        
        # Convert vectors to list
        self.ReferenceFileNames = getListFromVector(ReferenceFileNames)
        self.numberOfReferences = len(self.ReferenceFileNames)
        #directory with ProjMatchClasses
        self.ProjMatchDirs=[" "]
        self.LibraryDirs=[" "]
        self.DocFileInputAngles=[self.DocFileWithOriginalAngles]
        #ProjMatchRootName=[" "]
        
        for iterN in range(NumberofIterations):
            fnBaseIter = "%s/Iter_%02d/" % (self.WorkingDir, iterN + 1)
            self.ProjMatchDirs.append(fnBaseIter + self.ProjMatchDir)
            self.LibraryDirs.append( fnBaseIter + self.LibraryDir)
            self.DocFileInputAngles.append("%s%s.doc" % (fnBaseIter, self.docfile_with_current_angles))
        
        auxList = (self.numberOfReferences + 1) * [None]
        self.ProjectLibraryRootNames=[[None]]
        for iterN in range(NumberofIterations):
            fnBaseIter = "%s/Iter_%02d/" % (self.WorkingDir, iterN + 1)
            for refN in range(self.numberOfReferences):                
                auxList[refN + 1]= "%s%s_ref_%02d.stk" % (fnBaseIter, self.ProjectLibraryRootName, refN)
            self.ProjectLibraryRootNames.append(list(auxList))
                    
        self.ProjMatchRootNames=[[None]]
        for iterN in range(NumberofIterations):
            for refN in range(self.numberOfReferences):
                auxList[refN + 1]="%s/%s_ref_%02d.doc" % (self.ProjMatchDirs[iterN + 1], ProjMatchName, refN + 1)
            self.ProjMatchRootNames.append(list(auxList))
    
    
        #name of masked volumes
        #add dummy name so indexes start a 1
        self.maskedFileNamesIters = [[None]]
        for iterN in range(NumberofIterations):
            fnBaseIter = "%s/Iter_%02d/" % (self.WorkingDir, iterN + 1)
            for refN in range(self.numberOfReferences):
                auxList[refN + 1] = "%s%s_ref_%02d.vol" % (fnBaseIter, self.maskReferenceVolume, refN + 1)
            self.maskedFileNamesIters.append(list(auxList))
    
        ####################################################################
        #add initial reference, useful for many routines
        #NOTE THAT INDEXES START AT 1
        self.reconstructedFileNamesIters.append([None] + self.ReferenceFileNames)
        for iterN in range(NumberofIterations):
            fnBaseIter = "%s/Iter_%02d/" % (self.WorkingDir, iterN + 1)
            for refN in range(self.numberOfReferences):
                auxList[refN + 1] = "%s%s_ref_%02d.vol" % (fnBaseIter, self.ReconstructedVolume, refN + 1)
            self.reconstructedFileNamesIters.append(list(auxList))
    
        # Optimal angles from previous iteration or user-provided at the beginning
    #    global self.DocFileInputAngles
    #    self.DocFileInputAngles=[]
    #    aux=[self.DocFileWithOriginalAngles]*(self.numberOfReferences+1)#+1 fills the zero
    #    self.DocFileInputAngles.append([None] + aux)
    #    for iterN in range(NumberofIterations):
    #        for refN in range(self.numberOfReferences):
    #            auxList[refN + 1] = WorkingDir + "/Iter_" + \
    #                                      str(iterN + 1).zfill(2) + \
    #                                      '/' + \
    #                                      self.docfile_with_current_angles + \
    #                                      "_ref_" + str(refN + 1).zfill(2) + ".doc"
    #        self.DocFileInputAngles.append(list(auxList))
    
        self.docfile_with_current_anglesList=[None]
        for iterN in range(NumberofIterations):
            fnBaseIter = "%s/Iter_%02d/%s.doc" % (self.WorkingDir, iterN + 1, self.docfile_with_current_angles)
            self.docfile_with_current_anglesList.append(fnBaseIter)
    
        #parameter for projection matching
        self.Align2DIterNr          = [-1]+getListFromVector(Align2DIterNr,NumberofIterations)
        self.Align2dMaxChangeOffset = [-1]+getListFromVector(Align2dMaxChangeOffset,NumberofIterations)
        self.Align2dMaxChangeRot    = [-1]+getListFromVector(Align2dMaxChangeRot,NumberofIterations)
        self.AngSamplingRateDeg     = [-1]+getListFromVector(AngSamplingRateDeg,NumberofIterations)
        self.DiscardPercentage      = [-1]+getListFromVector(DiscardPercentage,NumberofIterations)
        self.DoAlign2D              = [False]+getBoolListFromVector(DoAlign2D,NumberofIterations)
        self.DoComputeResolution    = [False]+getBoolListFromVector(DoComputeResolution,NumberofIterations)
        self.DoSplitReferenceImages = [False]+getBoolListFromVector(DoSplitReferenceImages,NumberofIterations)
        self.InnerRadius            = [False]+getListFromVector(InnerRadius,NumberofIterations)
        self.MaxChangeInAngles      = [-1]+getListFromVector(MaxChangeInAngles,NumberofIterations)
        self.MaxChangeOffset        = [-1]+getListFromVector(MaxChangeOffset,NumberofIterations)
        self.MinimumCrossCorrelation= [-1]+getListFromVector(MinimumCrossCorrelation,NumberofIterations)
        self.OnlyWinner             = [False]+getBoolListFromVector(OnlyWinner,NumberofIterations)
        self.OuterRadius            = [False]+getListFromVector(OuterRadius,NumberofIterations)
        self.PerturbProjectionDirections = [False]+getBoolListFromVector(PerturbProjectionDirections,NumberofIterations)
        self.ReferenceIsCtfCorrected     = [-1]+getListFromVector(str(ReferenceIsCtfCorrected) + " True", NumberofIterations)
        self.ScaleNumberOfSteps          = [-1]+getListFromVector(ScaleNumberOfSteps,NumberofIterations)
        self.ScaleStep              = [-1]+getListFromVector(ScaleStep,NumberofIterations)
        self.Search5DShift          = [-1]+getListFromVector(Search5DShift,NumberofIterations)
        self.Search5DStep           = [-1]+getListFromVector(Search5DStep,NumberofIterations)
        self.SymmetryGroup          = [-1]+getListFromVector(SymmetryGroup,NumberofIterations)

    def otherActionsToBePerformedBeforeLoop(self):    
        #global OuterRadius, 
        global NumberOfCtfGroups
        
        if DoCtfCorrection:
            auxMD1 = MetaData(CTFDatName)
            auxMD2 = MetaData()
            auxMD2.aggregate(auxMD1, AGGR_COUNT,MDL_CTFMODEL,MDL_CTFMODEL,MDL_COUNT)
            NumberOfCtfGroups = auxMD2.size()
        else:
            NumberOfCtfGroups = 1
    
        _Parameters = {
              'DoDeleteWorkingDir':DoDeleteWorkingDir
            , 'ProjectDir':ProjectDir
            , 'WorkingDir':WorkingDir
            }
        command = 'deleteWorkingDirectory'
        _dataBase.insertCommand(command, _Parameters, 1)
    
        #Create directory
        _Parameters = {
              'Iter':0
            , 'ProjectDir':ProjectDir
            , 'WorkingDir':WorkingDir
            }
        command = 'createDir'
        _dataBase.insertCommand(command, _Parameters, 1)
    
        #Backup protocol file
        _Parameters = {
              'ProgName'  :sys.argv[0]
            , 'ProjectDir':ProjectDir
            , 'WorkingDir':WorkingDir
            }
        command = 'pm_make_backup_of_script_file'
        _dataBase.insertCommand(command, _Parameters, dataBase.dataBaseStruct.doAlways)#backup always
    
        #Check references and projections size match
        #Is already done in preconditions but I like to
        #run protocols from command line bypassing the gui
        _Parameters = {
              'ReferenceFileNames':self.ReferenceFileNames
            , 'SelFileName':SelFileName
            }
        command = 'checkVolumeProjSize'
        _dataBase.insertCommand(command, _Parameters, 1)
    
        #Check Option compatibility
        _Parameters = {
              'DoAlign2D':self.DoAlign2D[1]
            , 'DoCtfCorrection':DoCtfCorrection
            }
        command = 'checkOptionsCompatibility'
        _dataBase.insertCommand(command, _Parameters, 1)
    
    #    #Init Mask references radius
    #    _Parameters = {
    #          'OuterRadius':OuterRadius[1]
    #        , 'SelFileName':SelFileName
    #        }
    #    command = 'self.OuterRadius = initOuterRadius'
    #    _dataBase.insertCommand(command, _Parameters, 1)
    
        #7 make CTF groups
        _Parameters = {
                      'CTFDatName': CTFDatName
                    , 'CtfGroupDirectory': CtfGroupDirectory
                    , 'CtfGroupMaxDiff': CtfGroupMaxDiff
                    , 'CtfGroupMaxResol': CtfGroupMaxResol
                    , 'CtfGroupRootName': CtfGroupRootName
                    , 'DataArePhaseFlipped': DataArePhaseFlipped
                    , 'DoAutoCtfGroup': DoAutoCtfGroup
                    , 'DoCtfCorrection': DoCtfCorrection
                    , 'PaddingFactor': PaddingFactor
                    , 'SelFileName': SelFileName
                    , 'SplitDefocusDocFile': SplitDefocusDocFile
                    , 'WienerConstant': WienerConstant
                   }
        command = 'execute_ctf_groups'
        _VerifyFiles = []
        fnBase = os.path.join(self.CtfGroupDirectory, self.CtfGroupRootName)
        if DoCtfCorrection:            
            _VerifyFiles.append(fnBase+'Info.xmd')
            _VerifyFiles.append(fnBase+'_ctf.stk')
            _VerifyFiles.append(fnBase+'_wien.stk')
            _VerifyFiles.append(fnBase+'_split.doc')
        _VerifyFiles.append(fnBase+'_images.sel')
            
        _dataBase.insertCommand(command, _Parameters, 1,_VerifyFiles)
        #Create Initial angular file. Either fill it with zeros or copy input
        _Parameters = {
              'DocFileName':DocFileName
            , 'DocFileWithOriginalAngles':self.DocFileWithOriginalAngles
            , 'SelFileName':SelFileName
            }
        command = 'initAngularReferenceFile'
        _VerifyFiles = []
        _VerifyFiles.append(self.DocFileWithOriginalAngles)
        _dataBase.insertCommand(command, _Parameters, 1,_VerifyFiles)
    
        #Save all parameters in dict for future runs (this is done by database)
        #so far no parameter is being saved, but dummy=0
        _Parameters = {
          'SystemFlavour':SystemFlavour
        }
        command = 'self.saveParameters'
        _dataBase.insertCommand(command, _Parameters, 1)
        command = 'self.loadParameters'
        _dataBase.insertCommand(command, _Parameters, dataBase.dataBaseStruct.doAlways)
    
        #no entries will be save untill this commit
        print "commit databse"
        _dataBase.commit()
    
    def actionsToBePerformedInsideLoop(self, _log):
        
        for iterN in range(1, NumberofIterations + 1):
            #############conn.execute(sqlBegin + "MPI_ON" + sqlEnd)
            # create working dir
            _Parameters = {
                 'Iter':iterN
               , 'WorkingDir':WorkingDir
                }
            command = 'createDir'
            _dataBase.insertCommand(command, _Parameters, iterN)
    
            #Create directory with classes
            _Parameters = {
                  'path':self.ProjMatchDirs[iterN]
                }
            command = 'createDir2'
            _dataBase.insertCommand(command, _Parameters, 1)
        
            #Create directory with image libraries
            _Parameters = {
                  'path':self.LibraryDirs[iterN]
                }
            command = 'createDir2'
            _dataBase.insertCommand(command, _Parameters, 1)
            for refN in range(1, self.numberOfReferences + 1):
                ##############REMOVE SHUTIL.COPY
                # Mask reference volume
                _Parameters = {
                                      'DoMask'             : DoMask
                                    , 'DoSphericalMask'    : DoSphericalMask
                                    , 'maskedFileName'     : self.maskedFileNamesIters[iterN][refN]
                                    , 'maskRadius'         : MaskRadius
                                    , 'reconstructedFileName' : self.reconstructedFileNamesIters[iterN - 1][refN]
                                    , 'userSuppliedMask'   : MaskFileName
                                    }
                
                command = "execute_mask"
                _VerifyFiles = []
                _VerifyFiles.append(self.maskedFileNamesIters[iterN][refN])
                _VerifyFiles.append(self.maskedFileNamesIters[iterN][refN])
                _dataBase.insertCommand(command, _Parameters, iterN,_VerifyFiles)
    
                # angular_project_library
                _Parameters = {
                                      'AngSamplingRateDeg':self.AngSamplingRateDeg[iterN]
                                    , 'CtfGroupSubsetFileName':self.CtfGroupSubsetFileName
                                    , 'DoCtfCorrection': DoCtfCorrection
                                    , 'DocFileInputAngles':self.DocFileInputAngles[iterN-1]
                                    , 'DoParallel': DoParallel
                                    , 'DoRestricSearchbyTiltAngle':DoRestricSearchbyTiltAngle
                                    , 'MaxChangeInAngles':self.MaxChangeInAngles[iterN]
                                    , 'maskedFileNamesIter':self.maskedFileNamesIters[iterN][refN]
                                    , 'MpiJobSize':MpiJobSize
                                    , 'NumberOfMpiProcesses':NumberOfMpiProcesses
                                    , 'NumberOfThreads':NumberOfThreads
                                    , 'OnlyWinner':self.OnlyWinner[iterN]
                                    , 'PerturbProjectionDirections':self.PerturbProjectionDirections[iterN]
                                    , 'ProjectLibraryRootName':self.ProjectLibraryRootNames[iterN][refN]
                                    , 'SystemFlavour':SystemFlavour
                                    , 'SymmetryGroup':self.SymmetryGroup[iterN]
                                    , 'SymmetryGroupNeighbourhood':SymmetryGroupNeighbourhood
                                    , 'Tilt0':Tilt0
                                    , 'TiltF':TiltF
                                    }
    
                command = "angular_project_library"
                _VerifyFiles = []
                #file with projections
                auxFn=self.ProjectLibraryRootNames[iterN][refN]
                _VerifyFiles.append(auxFn)
                auxFn=auxFn[:-4]#remove extension
                #file with projection angles angles 
                _VerifyFiles.append(auxFn + ".doc")
                #file with sampling point neighbourhood 
                _VerifyFiles.append(auxFn + "_sampling.xmd")
                #file with sampling point neighbourhood for each ctf group, this is reduntant but useful
                
                for i in range (1,NumberOfCtfGroups+1):
                    _VerifyFiles.append(auxFn + "_group" + str(i).zfill(6) +"_sampling.xmd")
                            
                _dataBase.insertCommand(command, _Parameters, iterN,_VerifyFiles)
                # projectionMatching    
                _Parameters = {
                                      'AvailableMemory':AvailableMemory
                                    , 'CtfGroupRootName': CtfGroupRootName
                                    , 'CtfGroupDirectory': self.CtfGroupDirectory
                                    , 'DoComputeResolution':self.DoComputeResolution[iterN]
                                    , 'DoCtfCorrection': DoCtfCorrection
                                    , 'DoScale':DoScale
                                    , 'DoParallel': DoParallel
                                    , 'InnerRadius':self.InnerRadius[iterN]
                                    , 'MaxChangeOffset':self.MaxChangeOffset[iterN]
                                    , 'MpiJobSize':MpiJobSize
                                    , 'NumberOfCtfGroups':NumberOfCtfGroups
                                    , 'NumberOfMpiProcesses':NumberOfMpiProcesses
                                    , 'NumberOfThreads':NumberOfThreads
                                    , 'OuterRadius':self.OuterRadius[iterN]
                                    , 'PaddingFactor':PaddingFactor
                                    , 'ProjectLibraryRootName':self.ProjectLibraryRootNames[iterN][refN]
                                    , 'ProjMatchRootName':self.ProjMatchRootNames[iterN][refN]
                                    , 'ReferenceIsCtfCorrected':self.ReferenceIsCtfCorrected[iterN]
                                    , 'ScaleStep':self.ScaleStep[iterN]
                                    , 'ScaleNumberOfSteps':self.ScaleNumberOfSteps[iterN]
                                    , 'Search5DShift':self.Search5DShift[iterN]
                                    , 'Search5DStep':self.Search5DStep[iterN]
                                    , 'SystemFlavour':SystemFlavour
                                    }
    
                command = "projection_matching"
                _VerifyFiles = []
                #File with list of images and references
                _VerifyFiles.append(self.ProjMatchRootNames[iterN][refN] )
                for i in range (1,NumberOfCtfGroups+1):
                    _VerifyFiles.append(auxFn + "_group" + str(i).zfill(6) +"_sampling.xmd")
                _dataBase.insertCommand(command, _Parameters, iterN,_VerifyFiles)
                
    
            
            #assign the images to the different references based on the crosscorrelation coheficient
            #if only one reference it just copy the docfile generated in the previous step
            _Parameters = {
                           'DocFileInputAngles' : self.DocFileInputAngles[iterN]#Output file with angles
                         , 'NumberOfCtfGroups' : NumberOfCtfGroups
                         , 'ProjMatchRootName':self.ProjMatchRootNames[iterN]#LIST
                         , 'NumberOfReferences':self.numberOfReferences
                          }
            _VerifyFiles = []
            _VerifyFiles.append(self.DocFileInputAngles[iterN])
            command = "assign_images_to_references"
            _dataBase.insertCommand(command, _Parameters, iterN,_VerifyFiles)
    
            #align images, not possible for ctf groups
            for refN in range(1, self.numberOfReferences + 1):
                _Parameters = {
                           'Align2DIterNr':self.Align2DIterNr[iterN]#
                         , 'Align2dMaxChangeRot':self.Align2dMaxChangeRot[iterN]#
                         , 'Align2dMaxChangeOffset':self.Align2dMaxChangeOffset[iterN]#
                         , 'CtfGroupDirectory': self.CtfGroupDirectory#
                         , 'CtfGroupRootName': CtfGroupRootName#
                         , 'DiscardPercentage':self.DiscardPercentage[iterN]#
                         , 'DoAlign2D' : self.DoAlign2D[iterN]#
                         , 'DoComputeResolution':self.DoComputeResolution[iterN]
                         , 'DoCtfCorrection': DoCtfCorrection#
                         , 'DocFileInputAngles' : self.DocFileInputAngles[iterN]#
                         , 'DoParallel': DoParallel#
                         , 'DoSplitReferenceImages':self.DoSplitReferenceImages[iterN]#
                         , 'InnerRadius':self.InnerRadius[iterN]#
                         , 'MaxChangeOffset':self.MaxChangeOffset[iterN]#
                         , 'MinimumCrossCorrelation':self.MinimumCrossCorrelation[iterN]#
                         , 'NumberOfReferences':self.numberOfReferences#
                         , 'NumberOfCtfGroups' : NumberOfCtfGroups#
                         , 'NumberOfMpiProcesses':NumberOfMpiProcesses#
                         , 'NumberOfThreads':NumberOfThreads#
                         , 'PaddingFactor':PaddingFactor#
                         , 'ProjectLibraryRootName':self.ProjectLibraryRootNames[iterN][refN]#
                         , 'ProjMatchRootName':self.ProjMatchRootNames[iterN][refN]#
                         , 'refN':refN
                         , 'SystemFlavour':SystemFlavour#
                        }
    #, 'MpiJobSize':MpiJobSize#
    #, 'OuterRadius':OuterRadius[iterN]
                command = "angular_class_average"
                _VerifyFiles = []
                _VerifyFiles.append('ertertertert')
                _dataBase.insertCommand(command, _Parameters, iterN,_VerifyFiles)
                
                ##############REMOVE SHUTIL.COPY
                # Mask reference volume
                _Parameters = {
                                      'Align2dMaxChangeOffset':self.Align2dMaxChangeOffset[iterN]
                                    , 'Align2dMaxChangeRot':self.Align2dMaxChangeRot[iterN]
                                    , 'DoMask'             : DoMask
                                    , 'DoSphericalMask'    : DoSphericalMask
                                    , 'maskedFileName'     : self.maskedFileNamesIters[iterN][refN]
                                    , 'maskRadius'         : MaskRadius
                                    , 'reconstructedFileName' : self.reconstructedFileNamesIters[iterN - 1][refN]
                                    , 'userSuppliedMask'   : MaskFileName
                                    }
                
                command = "execute_mask"
                _VerifyFiles = []
                _VerifyFiles.append(self.maskedFileNamesIters[iterN][refN])
                _dataBase.insertCommand(command, _Parameters, iterN,_VerifyFiles)
    
    
    
    
    
    #class average
    #alig2d
    
    #reconstruct
    #resolution
    #
                
    ######################
    ######################
    ########################            
                #REMOVE
                #Create directory
            _Parameters = {
                 'Iter':iterN + 1
               , 'WorkingDir':WorkingDir
                }
            command = 'createDir'
            _dataBase.insertCommand(command, _Parameters, iterN)
            
            _Parameters = {'dummy':''}
            command = "shutil.copy('%s','%s');dummy" % (self.ReferenceFileNames[0], self.reconstructedFileNamesIters[iterN][refN])
            _VerifyFiles = []
            _VerifyFiles.append(self.reconstructedFileNamesIters[iterN][refN])
            _dataBase.insertCommand(command, _Parameters, iterN)
    
    ###            #delete DocFileInputAngles so I can use append style for metadata in DocFileInputAngles
    ###            _Parameters = {
    ###                           'FileName': DocFileInputAngles[iterN]
    ###                          ,'Verbose' : 1
    ###                           }
    ###            command = "deleteFile"
    ###            _dataBase.insertCommand(command, _Parameters, iterN)
    
    ###        command = "exit(1)"
    ###        _dataBase.insertCommand(command, _Parameters, iterN)
    
    
    ######################################
        _Parameters = {'dummy':''}
        command = "print 'ALL DONE';dummy"
        _dataBase.insertCommand(command, _Parameters, dataBase.dataBaseStruct.doAlways)
    
        _dataBase.commit()

