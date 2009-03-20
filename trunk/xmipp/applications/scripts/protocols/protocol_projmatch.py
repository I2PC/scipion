#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for projection matching
#
# Example use:
# ./xmipp_protocol_projmatch.py
#
# Authors: Roberto Marabini,
#          Sjors Scheres,    March 2008
#
#-----------------------------------------------------------------------------
# {section} Global parameters
#-----------------------------------------------------------------------------
# {file} Selfile with the input images:
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
SelFileName='all_images.sel'

# {file} {expert} Docfile with the input angles:
""" Do not provide anything if there are no angles yet. 
    In that case, the starting angles will be read from the image headers
    This docfile should be in newXmipp-style format (with filenames as comments)
    Note that all filenames in this docfile should be with absolute paths!
"""
DocFileName=''

# {file} Initial 3D reference map:
ReferenceFileName='reference.vol'

# Working subdirectory: 
""" This directory will be created if it doesn't exist, and will be used to store all output from this run. Don't use the same directory for multiple different runs, instead use a structure like run1, run2 etc. 
"""
WorkingDir='ProjMatch/run1'

# Delete working subdirectory if it already exists?
""" Just be careful with this option...
"""
DoDeleteWorkingDir=True

# Number of iterations to perform
NumberofIterations=10

# Resume at iteration
""" This option may be used to finish a previously performed run.
    Set to 1 to start a new run 
    Note: Do NOT delete working directory if this option is not set to 1
"""
ContinueAtIteration=1

# {expert} Save disc space by cleaning up intermediate files?
""" Be careful, many options of the visualization protocol will not work anymore, 
    since all class averages, selfiles etc will be deleted.
"""
CleanUpFiles=False

# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project. Often, each data set of a given sample has its own ProjectDir.
"""
ProjectDir='/home/scheres/work/test'

# {expert} Directory name for logfiles:
LogDir='Logs'

#-----------------------------------------------------------------------------
# {section} CTF correction
#-----------------------------------------------------------------------------

# Perform CTF correction?
""" If set to true, a CTF (amplitude and phase) corrected map will be refined,
    and the data will be processed in CTF groups.
    Note that you cannot combine CTF-correction with re-alignment of the classes.
"""
DoCtfCorrection=True

# {file} CTFDat file with CTF data:
""" The input selfile may be a subset of the images in the CTFDat file, but all 
    images in the input selfile must be present in the CTFDat file. This field is 
    obligatory if CTF correction is to be performed. 
    Note that this file should be positioned in the project directory, and that the
    image names and ctf parameter filenames should be in absolute paths.
"""
CTFDatName='all_images.ctfdat'

# Make CTF groups automatically?
""" Make CTF groups based on a maximum differences at a given resolution limit.
    If this option is set to false, a docfile with the defocus values where to 
    split the images in distinct defocus group has to be provided (see expert option below)
"""
DoAutoCtfGroup=True

# Maximum difference in CTF-values in one group
""" If the difference between the CTF-values up to the resolution limit specified 
    below is larger than the value given here, two images will be placed in 
    distinct CTF groups.
"""
CtfGroupMaxDiff=0.5

# Resolution limit (Ang) for CTF-differences in one group
""" Maximum resolution where to consider CTF-differences among different groups.
    One should use somewhat higher resolutions than those aimed for in the refinement.
"""
CtfGroupMaxResol=15

# {file} {expert} Docfile with defocus values where to split into groups
""" This field is obligatory if you do not want to make the CTF groups automatically.
    Note that the requested docfile can be made initially with the xmipp_ctf_group program,
    and then it can be edited manually to suit your needs. 
"""
SplitDefocusDocFile=''

# {expert} Padding factor
""" Application of CTFs to reference projections and of Wiener filter to class averages will be done using padded images.
    Use values larger than one to pad the images. Suggestion, use 1 for large image and 2 for small
"""
PaddingFactor=1.

# {expert} Wiener constant
""" Term that will be added to the denominator of the Wiener filter.
    In theory, this value is the inverse of the signal-to-noise ratio
    If a negative value is taken, the program will use a default value as in FREALIGN 
    (i.e. 10% of average sum terms over entire space) 
    see Grigorieff JSB 157 (2006) pp117-125
"""
WienerConstant=-1

# Images have been phase flipped?
DataArePhaseFlipped=True

# Is the initial reference map CTF (amplitude) corrected?
ReferenceIsCtfCorrected=True

#-----------------------------------------------------------------------------
# {section} Mask
#-----------------------------------------------------------------------------
# Mask reference volume?
""" Masking the reference volume will increase the signal to noise ratio.
    Do not provide a very tight mask.
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mask for details
"""
DoMask=True

# {expert} Show masked volume
""" Masked volume will be shown. 
    Do not set this option to true for non-interactive processing (jobs sent to queues)
"""
DisplayMask=False

# {file} Binary mask-file used to mask the reference volume
MaskFileName='mask.vol'

#-----------------------------------------------------------------------------
# {section} Projection Matching
#-----------------------------------------------------------------------------
# Perform projection Matching?
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Projection_matching and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_projection_matching for details
"""
DoProjectionMatching=True

# {expert} Show projection maching library and classes
""" Show average of projections. Do not set this option to true for non-interactive processing (jobs sent to queues)
"""
DisplayProjectionMatching=False

# Inner radius for rotational correlation:
""" In pixels from the image center
"""
InnerRadius=0

# Outer radius for rotational correlation
""" In pixels from the image center. Use a negative number to use the entire image.
    WARNING: this radius will be use for masking before computing resoution
"""
OuterRadius=64

# {expert} Available memory to store all references (Gb)
""" This is only for the storage of the references. If yuor memories so not fit in memory, the projection matching program will run MUCH slower. But, keep in mind that probably some additional memory is needed for the operating system etc.
"""
AvailableMemory='2'

# Angular sampling rate
"""Angular distance (in degrees) between neighboring projection  points
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "8 8 2 2 " 
    specifies 4 iterations, the first two set the value to 8 
    and the last two to 2. An alternative compact notation 
    is ("2x8 2x0", i.e.,
    2 iterations with value 8, and 2 with value 2).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
AngSamplingRateDeg='4x10 2x5 2x3 2x2'

# Angular search range 
"""Maximum change in rot & tilt  (in +/- degrees)
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "1000 1000 10 10 " 
    specifies 4 iterations, the first two set the value to 1000 (no restriction)
    and the last two to 10degrees. An alternative compact notation 
    is ("2x1000 2x10", i.e.,
    2 iterations with value 1000, and 2 with value 10).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
MaxChangeInAngles='4x1000 2x20 2x9 2x6'

# {expert} Perturb projection directions?
""" If set to true, this option will result to a Gaussian perturbation to the 
    evenly sampled projection directions of the reference library. 
    This may serve to decrease the effects of model bias.
"""
PerturbProjectionDirections=False

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
MaxChangeOffset='1000 '

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
Search5DShift='4x5 0'

# {expert} Step size for 5D translational search
""" Provide a sequence of numbers (for instance, "2 2 1 1" specifies 4 iterations,
    the first two set the value to 2, then two with 1 pixel.
    An alternative compact notation is ("2x2 2x1", i.e.,
    2 iterations with value 2, and 2 with value 1).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
    
"""
Search5DStep='2'

# {expert} Restrict tilt angle search?
DoRetricSearchbyTiltAngle=False

# {expert} Lower-value for restricted tilt angle search
Tilt0=40

# {expert} Higher-value for restricted tilt angle search
TiltF=90

# Symmetry group
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry
    for a description of the symmetry groups format
    If no symmetry is present, give c1
"""
SymmetryGroup='c1'

# Discard images with ccf below
""" Provide a sequence of numbers (for instance, "0.3 0.3 0.5 0.5" specifies 4 iterations,
    the first two set the value to 0.3, then two with 0.5.
    An alternative compact notation would be ("2x0.3 2x0.5").
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
    Set to -1 to prevent discarding any images
"""    
MinimumCrossCorrelation='-1'

# Discard percentage of images with ccf below
""" Provide a sequence of numbers (for instance, "20 20 10 10" specifies 4 iterations,
    the first two set the value to 20%, then two with 10%
    An alternative compact notation would be ("2x20 2x10").
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
    Set to zero to prevent discarding any images
"""    
DiscardPercentage='10'

# {expert} Additional options for Projection_Matching
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Projection_matching and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_projection_matching for details
    try -Ri xx -Ro yy for restricting angular search (xx and yy are
    the particle inner and outter radius)
    
"""
ProjMatchingExtra=''

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
DoAlign2D=' 0'

# {expert} Number of align2d iterations:
""" Use at least 3 iterations
"""
Align2DIterNr=4

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
Align2dMaxChangeOffset='2x1000 2x10'

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
Align2dMaxChangeRot='2x1000 2x20'

#-----------------------------------------------------------------------------
# {section} 3D Reconstruction
#-----------------------------------------------------------------------------
# Perform 3D Reconstruction?
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Wbp and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_wbp and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Art
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Fourier
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_fourier
        for details
"""
DoReconstruction=True

# {expert} Display reconstructed volume?
DisplayReconstruction=False

# {list}|fourier|art|wbp| Reconstruction method
""" Choose between wbp, art or fourier
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "wbp wbp wbp art " 
    specifies 4 iterations, the first three set the value to wbp (no restriction)
    and the last  to art. An alternative compact notation 
    is ("3xwbp 1xart", i.e.,
    3 iterations with wbp, and 1 with art).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
ReconstructionMethod='fourier'

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
ARTLambda='0.2'

# {expert} Additional reconstruction parameters for ART
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Art
        for details
"""
ARTReconstructionExtraCommand='-k 0.5 -n 10 '

# Initial maximum frequency used by reconstruct fourier
""" This number os only used in the first iteration. 
    From then on, it will be set to resolution computed in the resolution section
"""
FourierMaxFrequencyOfInterest='0.25'

# {expert} Additional reconstruction parameters for WBP
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Wbp and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_wbp and
        for details
"""
WBPReconstructionExtraCommand=' '

# {expert} Additional reconstruction parameters for Fourier
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Fourier and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_Fourier and
        for details
	-thr_width 
"""
FourierReconstructionExtraCommand=' '

#-----------------------------------------------------------------------------
# {section} Compute Resolution
#-----------------------------------------------------------------------------
# Compute resolution?
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Resolution fo details
"""
DoComputeResolution=True

# Split references averages
"""In theory each reference average should be splited
   in two when computing the resolution. In this way each
   projection direction will be represented in each of the
   subvolumes used to compute the resolution. A much faster
   but less accurate approach is to split the 
   proyection directions in two but not the averages. We
   recomend the first approach for small volumes and the second for
   large volumes (expecially when using small angular
   sampling rates.
   IMPORTANT: the second option has ONLY been implemented for FOURIER
   reconstruction method. Other reconstruction methods require this
   flag to be set to True
"""
DoSplitReferenceImages=True


# Pixel size (in Ang.)
""" This will make that the X-axis in the resolution plots has units 1/Angstrom
"""
ResolSam=2.8

# {expert} Display resolution?
DisplayResolution=False

#-----------------------------------------------------------------------------
# {section} Low-pass filtering
#-----------------------------------------------------------------------------
# Provide your own filtration frequecy?
"""By default the volume will be filtered at a frecuency equal to
   the  resolution computed with resolution_fsc
   plus a constant provided by the user in the next
   input box. If this option is set to true then the
   filtration will be made at the constant value provided by
   the user WITHOUT adding the resolution computed by resolution_fsc.
"""
SetResolutiontoZero=False


# Constant to by add to the estimate resolution
""" The meaning of this field depends on the previous flag.
    If set as False then  the volume will be filtered at a frecuency equal to
    the  resolution computed with resolution_fsc
    plus the value provided in this field (units pixel^-1)
    If set to False  the volume will be filtered at the resolution
    provided in this field.
    Set this value to any value larger than .5 to avoid filtration.
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, ".15 .15 .1 .1" 
    specifies 4 iterations, the first two set the constant to .15
    and the last two to 0.1. An alternative compact notation 
    is ("2x.15 2x0.1", i.e.,
    4 iterations with value 0.15, and three with value .1).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
ConstantToAddToFiltration='0.1'

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Number of (shared-memory) threads?
""" This option provides shared-memory parallelization on multi-core machines. 
    It does not require any additional software, other than xmipp
"""
NumberOfThreads=1

# distributed-memory parallelization (MPI)?
""" This option provides distributed-memory parallelization on multi-node machines. 
    It requires the installation of some MPI flavour, possibly together with a queueing system
"""
DoParallel=True

# Number of MPI processes to use:
NumberOfMpiProcesses=5

# minumum size of jobs in mpi processe. Set to 1 for large images (e.g. 500x500) and to 10 for small images (e.g. 100x100)
MpiJobSize='10'

# MPI system Flavour 
""" Depending on your queuing system and your mpi implementation, different mpirun-like commands have to be given.
    Ask the person who installed your xmipp version, which option to use. Or read: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/ParallelPage. The following values are available: 
   SLURM-MPICH       : BSCs MareNostrum, LaPalma etc
   TORQUE-OPENMPI    : Crunchy
   SGE-OPENMPI       : Cluster at imp.ac.at
   PBS               : Vermeer and FinisTerrae
   XMIPP_MACHINEFILE : Environment variable $XMIPP_MACHINEFILE points to 
                       machinefile
   HOME_MACHINEFILE  : machinefile is called $HOME/machines.dat
   Leave it black    : Run locally, most personal computers
"""
SystemFlavour=''

#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript='visualize_projmatch.py'

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#-----------------------------------------------------------------------------
#Do not change these variables
ReferenceVolumeName='reference_volume.vol'
LibraryDir = "ReferenceLibrary"
ProjectLibraryRootName= LibraryDir + "/ref"
ProjMatchDir = "ProjMatchClasses"
ProjMatchRootName= ProjMatchDir + "/proj_match"
ForReconstructionSel="reconstruction.sel"
ForReconstructionDoc="reconstruction.doc"
MultiAlign2dSel="multi_align2d.sel"
DocFileWithOriginalAngles='original_angles.doc'
docfile_with_current_angles='current_angles.doc'
FilteredReconstruction="filtered_reconstruction"
ReconstructedVolume="reconstruction"
OutputFsc="resolution.fsc"
CtfGroupDirectory="CtfGroups"
CtfGroupRootName="ctf"
CtfGroupSubsetFileName="ctf_groups_subset_docfiles.sel"

class projection_matching_class:

   #init variables
   
   def __init__(self,
                _NumberofIterations,
                _ContinueAtIteration,
                _CleanUpFiles,
                _DoMask, 
                _DisplayMask,
                _ReferenceFileName,
                _MaskFileName,
                _DoProjectionMatching,
                _DisplayProjectionMatching,
                _AngSamplingRateDeg,
                _PerturbProjectionDirections,
                _DoRetricSearchbyTiltAngle,
                _Tilt0,
                _TiltF,
                _ProjMatchingExtra,
                _MaxChangeOffset,
                _MaxChangeInAngles,
                _MinimumCrossCorrelation,
                _DiscardPercentage,
                _DoAlign2D,
                _InnerRadius,
                _OuterRadius,
                _Search5DShift,
                _Search5DStep,
                _AvailableMemory,
                _Align2DIterNr,
                _Align2dMaxChangeOffset,
                _Align2dMaxChangeRot,
                _DisplayReconstruction,
                _DisplayResolution,
                _DoReconstruction,
                _ReconstructionMethod,
                _ARTLambda,
                _ARTReconstructionExtraCommand,
                _WBPReconstructionExtraCommand,
                _FourierReconstructionExtraCommand,
                _FourierMaxFrequencyOfInterest,
                _DoComputeResolution,
                _DoSplitReferenceImages,
                _ResolSam,
                _SelFileName,
                _DocFileName,
                _DoCtfCorrection,
                _CTFDatName,
                _WienerConstant,
                _DoAutoCtfGroup,
                _CtfGroupMaxDiff,
                _CtfGroupMaxResol,
                _SplitDefocusDocFile,
                _PaddingFactor,
                _DataArePhaseFlipped,
                _ReferenceIsCtfCorrected,
                _WorkingDir,
                _ProjectDir,
                _LogDir,
                _DoParallel,
                _MyNumberOfMpiProcesses,
                _MySystemFlavour,
                _MyMpiJobSize,
                _MyNumberOfThreads,
                _SymmetryGroup,
                _SetResolutiontoZero,
                _ConstantToAddToFiltration
                ):

       # Import libraries and add Xmipp libs to default search path
       import os,sys,shutil
       scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
       sys.path.append(scriptdir)
       import arg,log,logging,selfile
       import launch_job

       self._CleanUpFiles=_CleanUpFiles
       self._WorkingDir=os.getcwd()+'/'+_WorkingDir
       self._SelFileName=_SelFileName
       selfile_without_ext=(os.path.splitext(str(os.path.basename(self._SelFileName))))[0]
       self._ReferenceFileName=os.path.abspath(_ReferenceFileName)
       self._MaskFileName=os.path.abspath(_MaskFileName)
       self._DoMask=_DoMask
       self._DoProjectionMatching=_DoProjectionMatching
       self._DisplayProjectionMatching=_DisplayProjectionMatching
       self._DoRetricSearchbyTiltAngle=_DoRetricSearchbyTiltAngle
       self._PerturbProjectionDirections=_PerturbProjectionDirections
       self._Tilt0=_Tilt0
       self._TiltF=_TiltF
       self._ProjMatchingExtra=_ProjMatchingExtra
       self._DisplayMask=_DisplayMask
       self._ProjectDir=_ProjectDir
       self._InnerRadius=_InnerRadius
       self._AvailableMemory=_AvailableMemory
       self._Align2DIterNr=_Align2DIterNr
       self._DisplayReconstruction=_DisplayReconstruction
       self._DisplayResolution=_DisplayResolution
       self._DoReconstruction=_DoReconstruction
       self._DoComputeResolution=_DoComputeResolution
       self._DoSplitReferenceImages=_DoSplitReferenceImages
       self._ResolSam=_ResolSam
       self._DoCtfCorrection=_DoCtfCorrection
       self._WienerConstant=_WienerConstant
       self._DoAutoCtfGroup=_DoAutoCtfGroup
       self._CtfGroupMaxDiff=_CtfGroupMaxDiff
       self._CtfGroupMaxResol=_CtfGroupMaxResol
       self._SplitDefocusDocFile =''  
       if(len(_SplitDefocusDocFile) > 1):
            self._SplitDefocusDocFile=os.path.abspath(_SplitDefocusDocFile)
       self._DocFileName =''  
       if(len(_DocFileName) > 1):
            self._DocFileName=os.path.abspath(_DocFileName)
       self._PaddingFactor=PaddingFactor
       self._DataArePhaseFlipped=_DataArePhaseFlipped
       self._DoParallel=_DoParallel
       self._MyNumberOfMpiProcesses=_MyNumberOfMpiProcesses
       self._SymmetryGroup=_SymmetryGroup
       self._ARTReconstructionExtraCommand=_ARTReconstructionExtraCommand
       self._WBPReconstructionExtraCommand=_WBPReconstructionExtraCommand
       self._FourierReconstructionExtraCommand=_FourierReconstructionExtraCommand
       self._SetResolutiontoZero=_SetResolutiontoZero
       # if we are not starting at the first iteration
       # globalFourierMaxFrequencyOfInterest must be computed
       # untill I fix this properlly let us set it at max_frequency,
      
       if _ContinueAtIteration==1:
           globalFourierMaxFrequencyOfInterest=float(_FourierMaxFrequencyOfInterest)
       else:
           globalFourierMaxFrequencyOfInterest=0.5
       self._MySystemFlavour=_MySystemFlavour
       
       self._MyMpiJobSize =_MyMpiJobSize
       self._MyNumberOfThreads =_MyNumberOfThreads
       self._user_suplied_ReferenceVolume=self._ReferenceFileName

       # Set up logging
       self._mylog=log.init_log_system(_ProjectDir,
                                       _LogDir,
                                       sys.argv[0],
                                       _WorkingDir)
                                      
       # Uncomment next line to get Debug level logging
       self._mylog.setLevel(logging.DEBUG)
       self._mylog.debug("Debug level logging enabled")
                                      
       _NumberofIterations +=1;
       if _ContinueAtIteration!=1 and DoDeleteWorkingDir==True:
          print "You can not delete the working directory"
          print " and start at iteration", _ContinueAtIteration
          exit(1)
       if (DoDeleteWorkingDir): 
          delete_working_directory(self._mylog,self._WorkingDir)
       else:
          self._mylog.info("Skipped DoDeleteWorkingDir") 
       if ReconstructionMethod!='fourier' and not _DoSplitReferenceImages:
          print "DoSplitReferenceImages must be set to True is"
          print " reconstructin method is not fourier"
          exit(1)
       create_working_directory(self._mylog,self._WorkingDir)
       log.make_backup_of_script_file(sys.argv[0],self._WorkingDir)
       
       # Create a selfile with absolute pathname in the WorkingDir
       mysel=selfile.selfile()
       mysel.read(_SelFileName)
       newsel=mysel.make_abspath()
       self._SelFileName=os.path.abspath(self._WorkingDir + '/' + _SelFileName)
       newsel.write(self._SelFileName)
       
       # For ctf groups, also create a CTFdat file with absolute pathname in the WorkingDir
       if (self._DoCtfCorrection):
          import ctfdat
          myctfdat=ctfdat.ctfdat()
          myctfdat.read(_CTFDatName)
          newctfdat=myctfdat.make_abspath()
          self._CTFDatName=os.path.abspath(self._WorkingDir + '/' + _CTFDatName)
          newctfdat.write(self._CTFDatName)

       # Set self._OuterRadius
       if (_OuterRadius < 0):
          xdim,ydim=newsel.imgSize()
          self._OuterRadius = (xdim/2) - 1 
          comment = " Outer radius set to: " + str(self._OuterRadius)
          print '* ' + comment
          self._mylog.info(comment)
       else:   
          self._OuterRadius=_OuterRadius

       # Create a docfile with the current angles in the WorkingDir
       if (self._DocFileName==''):

          params=' -i ' + self._SelFileName + \
                 ' -o ' + self._WorkingDir + '/' + \
                 DocFileWithOriginalAngles
          launch_job.launch_job("xmipp_header_extract",
                                params,
                                self._mylog,
                                False,1,1,'')
       else:
          command = "copy" , self._DocFileName ,  self._WorkingDir + '/' + DocFileWithOriginalAngles
          self._mylog.info(command)
          shutil.copy(self._DocFileName, self._WorkingDir + '/' + DocFileWithOriginalAngles)

       # Change to working dir
       os.chdir(self._WorkingDir)
       self._SelFileName=self._WorkingDir+'/'+\
                         str(os.path.basename(self._SelFileName))

       # Make CTF groups
       if (self._DoCtfCorrection):
          self._NumberOfCtfGroups=execute_ctf_groups(self._mylog,
                                                     self._SelFileName,
                                                     self._CTFDatName,
                                                     self._PaddingFactor,
                                                     self._DataArePhaseFlipped,
                                                     self._WienerConstant,
                                                     self._DoAutoCtfGroup,
                                                     self._CtfGroupMaxDiff,
                                                     self._CtfGroupMaxResol,
                                                     self._SplitDefocusDocFile)
       else:
          self._NumberOfCtfGroups=1

       ##
       ##LOOP
       ##
       #output of reconstruction cycle
       #first value given by user
       #these names are the input of the mask program
       #in general is the output of the reconstruction plus filtration
       self._ReconstructedVolume=[]
       fill_name_vector("",
                        self._ReconstructedVolume,
                        _NumberofIterations,
                        ReconstructedVolume)
                        
       self._ReconstructedandfilteredVolume=[]
       fill_name_vector(self._user_suplied_ReferenceVolume,
                        self._ReconstructedandfilteredVolume,
                        _NumberofIterations,
                        FilteredReconstruction)

       # Optimal angles from previous iteration or user-provided at the beginning
       self._DocFileInputAngles=[]
       fill_name_vector('../'+DocFileWithOriginalAngles,
                        self._DocFileInputAngles,
                        _NumberofIterations+1,
                        docfile_with_current_angles)

       # Reconstructed and filtered volume of n-1 after masking called reference volume
       self._ReferenceVolume=[]
       fill_name_vector("",
                        self._ReferenceVolume,
                        _NumberofIterations,
                        ReferenceVolumeName)

       for _iteration_number in range(_ContinueAtIteration, _NumberofIterations):
          debug_string =  "ITERATION: " +  str(_iteration_number)
          print "*", debug_string
          self._mylog.info(debug_string)

          # Never allow DoAlign2D and DoCtfCorrection together
          if (arg.getComponentFromVector(_DoAlign2D,_iteration_number-1)==1 and
              self._DoCtfCorrection):
             error_message="You cannot realign classes AND perform CTF-correction. Switch either of them off!"
             self._mylog.error(error_message)
             print error_message
             exit(1)

          # Create working dir for this iteration and go there
          Iteration_Working_Directory=self._WorkingDir+'/Iter_'+\
                                      str(_iteration_number)
          create_working_directory(self._mylog,Iteration_Working_Directory)
          os.chdir(Iteration_Working_Directory)

          # Mask reference volume
          execute_mask(_DoMask,
                       self._mylog,
                       self._ProjectDir,
                       self._ReconstructedandfilteredVolume[_iteration_number],#in
                       self._MaskFileName,
                       self._DisplayMask,
                       _iteration_number,
                       self._ReferenceVolume[_iteration_number])#out

          if (_DoProjectionMatching):
             # Parameters for projection matching
             self._AngSamplingRateDeg=arg.getComponentFromVector(_AngSamplingRateDeg,\
                                                           _iteration_number-1)
             self._MaxChangeOffset=arg.getComponentFromVector(_MaxChangeOffset,\
                                                           _iteration_number-1)
             self._MaxChangeInAngles=arg.getComponentFromVector(_MaxChangeInAngles,\
                                                           _iteration_number-1)
             self._Search5DShift=arg.getComponentFromVector(_Search5DShift,\
                                                           _iteration_number-1)
             self._Search5DStep=arg.getComponentFromVector(_Search5DStep,\
                                                           _iteration_number-1)
             self._MinimumCrossCorrelation=arg.getComponentFromVector(_MinimumCrossCorrelation,\
                                                           _iteration_number-1)
             self._DiscardPercentage=arg.getComponentFromVector(_DiscardPercentage,\
                                                           _iteration_number-1)
             self._DoAlign2D=arg.getComponentFromVector(_DoAlign2D,\
                                                           _iteration_number-1)
             self._Align2dMaxChangeOffset=arg.getComponentFromVector(_Align2dMaxChangeOffset,\
                                                           _iteration_number-1)
             self._Align2dMaxChangeRot=arg.getComponentFromVector(_Align2dMaxChangeRot,\
                                                           _iteration_number-1)

             # Initial reference is CTF-amplitude corrected?
             if ( (_iteration_number == 1) and (_ReferenceIsCtfCorrected==False) ):
                self._ReferenceIsCtfCorrected=False
             else: 
                self._ReferenceIsCtfCorrected=True

             execute_projection_matching(self._mylog,
                                         self._ProjectDir,
                                         self._ReferenceVolume[_iteration_number],
                                         self._MaskFileName,
                                         self._DocFileInputAngles[_iteration_number],
                                         self._DocFileInputAngles[_iteration_number+1],
                                         self._DoCtfCorrection,
                                         self._NumberOfCtfGroups,
                                         self._WienerConstant,
                                         self._PaddingFactor,
                                         self._ReferenceIsCtfCorrected,
                                         self._AngSamplingRateDeg,
                                         self._PerturbProjectionDirections,
                                         self._DoRetricSearchbyTiltAngle,
                                         self._Tilt0,
                                         self._TiltF,
                                         self._InnerRadius,
                                         self._OuterRadius,
                                         self._Search5DShift,
                                         self._Search5DStep,
                                         self._MaxChangeOffset, 
                                         self._MaxChangeInAngles,
                                         self._ProjMatchingExtra,
                                         self._MinimumCrossCorrelation,
                                         self._DiscardPercentage,
                                         self._DisplayProjectionMatching,
                                         self._DoParallel,
                                         self._MyNumberOfMpiProcesses,
                                         self._MyNumberOfThreads,
                                         self._MySystemFlavour,
                                         self._MyMpiJobSize,
                                         self._WorkingDir,
                                         self._SymmetryGroup,
                                         self._AvailableMemory,
                                         self._DoComputeResolution,
                                         self._DoSplitReferenceImages,
                                         self._DoAlign2D,
                                         self._Align2DIterNr,
                                         self._Align2dMaxChangeOffset,
                                         self._Align2dMaxChangeRot,
                                         _iteration_number
                                         )
          else:
             self._mylog.info("Skipped ProjectionMatching") 


          # Make a new selfile excluding the images that were possibly discarded by the user
          command='cat ' + MultiAlign2dSel + \
                        ' | grep -v ' +  ProjectLibraryRootName + \
                        ' | grep -v ref.xmp ' + \
                        ' | grep -v \ -1 >' + ForReconstructionSel
          self._mylog.info(command)
          os.system(command)

          self._ReconstructionMethod=arg.getComponentFromVector(_ReconstructionMethod,\
                                                        _iteration_number-1)
          self._ARTLambda=arg.getComponentFromVector(_ARTLambda,\
                                                        _iteration_number-1)
          if (_DoReconstruction):
             execute_reconstruction(self._mylog, 
                                    self._ARTReconstructionExtraCommand,
                                    self._WBPReconstructionExtraCommand,
                                    self._FourierReconstructionExtraCommand,
                                    _iteration_number,
                                    self._DisplayReconstruction,
                                    self._DoParallel,
                                    self._MyNumberOfMpiProcesses,
                                    self._MyNumberOfThreads,
                                    self._MySystemFlavour,
                                    self._MyMpiJobSize,
                                    self._ReconstructionMethod,
                                    globalFourierMaxFrequencyOfInterest,
                                    self._ARTLambda,
                                    self._SymmetryGroup,
                                    self._ReconstructedVolume[_iteration_number],
                                    self._DoComputeResolution,
                                    self._DoSplitReferenceImages
                                    )
          else:
             self._mylog.info("Skipped Reconstruction") 
          
          if (_DoComputeResolution):
              filter_frequence=execute_resolution(self._mylog,
                                                  self._ARTReconstructionExtraCommand,
                                                  self._WBPReconstructionExtraCommand,
                                                  self._FourierReconstructionExtraCommand,
                                                  self._ReconstructionMethod,
                                                  globalFourierMaxFrequencyOfInterest,
                                                  _iteration_number,
                                                  self._DisplayReconstruction,
                                                  self._ResolSam,
                                                  self._DoParallel,
                                                  self._MyNumberOfMpiProcesses,
                                                  self._MyNumberOfThreads,
                                                  self._MySystemFlavour,
                                                  self._MyMpiJobSize,
                                                  self._SymmetryGroup,
                                                  self._DisplayResolution,
                                                  self._ReconstructedVolume[_iteration_number],
                                                  self._ARTLambda,
                                                  self._OuterRadius
                                                  )
          else:
	     filter_frequence=0
             self._mylog.info("Skipped Resolution") 
          
          self._ConstantToAddToFiltration=arg.getComponentFromVector(\
                                               ConstantToAddToFiltration,\
                                                  _iteration_number-1)

          globalFourierMaxFrequencyOfInterest=filter_at_given_resolution(_DoComputeResolution,
                                     self._mylog, 
                                     _iteration_number,
                                     self._SetResolutiontoZero,
                                     self._ConstantToAddToFiltration,
                                     filter_frequence,
                                     self._ReconstructedVolume[_iteration_number],
                                     self._ReconstructedandfilteredVolume[1+_iteration_number],
                                     self._MySystemFlavour
                                     )

          # Remove all class averages and reference projections
          if (self._CleanUpFiles):
             execute_cleanup(self._mylog,
                             True,
                             True)


#------------------------------------------------------------------------
#delete_working directory
#------------------------------------------------------------------------
def delete_working_directory(_mylog,_WorkingDir):
    import os
    import shutil
    print '*********************************************************************'
    print '* Delete working directory tree'
    _mylog.info("Delete working directory tree")

    if os.path.exists(_WorkingDir):
       shutil.rmtree(_WorkingDir)
       
#------------------------------------------------------------------------
#create_working directory
#------------------------------------------------------------------------
def create_working_directory(_mylog,_WorkingDir):
    import os
    print '*********************************************************************'
    print '* Create directory ' + _WorkingDir 
    _mylog.info("Create working directory " + _WorkingDir )

    if not os.path.exists(_WorkingDir):
       os.makedirs(_WorkingDir)
    # Also create subdirectories
    if not os.path.exists(_WorkingDir + "/" + LibraryDir):
       os.makedirs(_WorkingDir + "/" + LibraryDir)
    if not os.path.exists(_WorkingDir + "/" + ProjMatchDir):
       os.makedirs(_WorkingDir + "/" + ProjMatchDir)

#------------------------------------------------------------------------
#make ctf groups
#------------------------------------------------------------------------
def execute_ctf_groups (_mylog,
                        _InPutSelfile,
                        _CtfDatFile,
                        _PaddingFactor,
                        _DataArePhaseFlipped,
                        _WienerConstant,
                        _DoAutoCtfGroup,
                        _CtfGroupMaxDiff,
                        _CtfGroupMaxResol,
                        _SplitDefocusDocFile):

   import os,glob,sys
   import utils_xmipp
   import launch_job

   if not os.path.exists(CtfGroupDirectory):
      os.makedirs(CtfGroupDirectory)

   print '*********************************************************************'
   print '* Make CTF groups'
   command=' -i '    + _InPutSelfile + \
           ' -ctfdat ' + _CtfDatFile + \
           ' -o ' + CtfGroupDirectory + '/' + CtfGroupRootName + \
           ' -wiener -wc ' + str(_WienerConstant) + \
           ' -pad ' + str(_PaddingFactor)
   if (_DataArePhaseFlipped):
      command += ' -phase_flipped '
   if (_DoAutoCtfGroup):
      command += ' -error ' + str(_CtfGroupMaxDiff) + \
                 ' -resol ' + str(_CtfGroupMaxResol)
   else:
      if (len(_SplitDefocusDocFile) > 0):
         command += ' -split ' + _SplitDefocusDocFile
      else:
         message = "Error: for non-automated ctf grouping, please provide a docfile!"
         print '* ',message
         _mylog.info(message)
         sys.exit()

   launch_job.launch_job("xmipp_ctf_group",
                         command,
                         _mylog,
                         False,1,1,'')

   wildcardname=utils_xmipp.composeWildcardFileName(CtfGroupDirectory + '/' + CtfGroupRootName+'_group','ctf')
   ctflist=glob.glob(wildcardname)
   return len(ctflist)

#------------------------------------------------------------------------
#execute_mask
#------------------------------------------------------------------------
def execute_mask(_DoMask,
                 _mylog,
                 _ProjectDir,
                 _ReferenceFileName,
                 _MaskFileName,
                 _DisplayMask,
                 _iteration_number,
                 _ReferenceVolume):
   import os,shutil
   import launch_job
   _mylog.debug("execute_mask")
   if(_iteration_number==1):
      InPutVolume=_ReferenceFileName
   else:   
      InPutVolume=_ReferenceFileName+".vol"
   if (_DoMask):
       MaskVolume =_MaskFileName
       MaskedVolume=_ReferenceVolume
       print '*********************************************************************'
       print '* Mask the reference volume'
       command=' -i '    + InPutVolume + \
               ' -o '    + _ReferenceVolume + \
               ' -mask ' + MaskVolume 
       launch_job.launch_job("xmipp_mask",
                             command,
                             _mylog,
                             False,1,1,'')

       if _DisplayMask==True:
          command='xmipp_show -vol '+ MaskedVolume +' -w 10 &'
          print '*********************************************************************'
          print '* ',command
          _mylog.info(command)
          os.system(command)
          
   else:
       shutil.copy(InPutVolume,_ReferenceVolume)
       _mylog.info("Skipped Mask")
       _mylog.info("cp" + InPutVolume +\
                   " "  + _ReferenceVolume )
       print '*********************************************************************'
       print '* Skipped Mask'

#------------------------------------------------------------------------
#execute_projection_matching
#------------------------------------------------------------------------
def execute_projection_matching(_mylog,
                                _ProjectDir,
                                _ReferenceVolume,
                                _MaskFileName,
                                _InputDocFileName,
                                _OutputDocFileName,
                                _DoCtfCorrection,
                                _NumberOfCtfGroups,
                                _WienerConstant,
                                _PaddingFactor,
                                _ReferenceIsCtfCorrected,
                                _AngSamplingRateDeg,
                                _PerturbProjectionDirections,
                                _DoRetricSearchbyTiltAngle,
                                _Tilt0,
                                _TiltF,
                                _Ri,
                                _Ro,
                                _Search5DShift,
                                _Search5DStep,
                                _MaxChangeOffset,
                                _MaxChangeInAngles,
                                _ProjMatchingExtra,
                                _MinimumCrossCorrelation,
                                _DiscardPercentage,
                                _DisplayProjectionMatching,
                                _DoParallel,
                                _MyNumberOfMpiProcesses,
                                _MyNumberOfThreads,
                                _MySystemFlavour,
                                _MyMpiJobSize,
                                _WorkingDir,
                                _SymmetryGroup,
                                _AvailableMemory,
                                _DoComputeResolution,
                                _DoSplitReferenceImages,
                                _DoAlign2D,
                                _Align2DIterNr,
                                _Align2dMaxChangeOffset,
                                _Align2dMaxChangeRot,
                                _iteration_number):
                                           
   _mylog.debug("execute_projection_matching")
   import os, shutil, string, glob, math
   import launch_job, selfile, docfiles, utils_xmipp

   if (_DoCtfCorrection):
      # To use -add_to in angular_class_average correctly, 
      # make sure there are no proj_match_class* files from previous runs. 
      print ' * CleanUp: deleting directory '+ ProjMatchDir
      os.system(' rm -r ' + ProjMatchDir)
      os.makedirs(ProjMatchDir)
      # Create docfiles for each defocus group and corresponding selfile containing all of them      
      make_subset_docfiles(_mylog,
                           _InputDocFileName,
                           _NumberOfCtfGroups)
      
   # Project all references
   print '*********************************************************************'
   print '* Create projection library'
   parameters=' -i '                   + _ReferenceVolume + \
              ' -experimental_images ' +  _InputDocFileName + \
              ' -o '                   + ProjectLibraryRootName + \
              ' -sampling_rate '       + _AngSamplingRateDeg  + \
              ' -sym '                 + _SymmetryGroup + 'h' + \
              ' -compute_neighbors '

   if ( string.atof(_MaxChangeInAngles) < 181.):
      parameters+= \
              ' -angular_distance '    + str(_MaxChangeInAngles)
   else:
      parameters+= \
              ' -angular_distance -1'

   if (_PerturbProjectionDirections):
      # Just follow Roberto's suggestion
      perturb=math.sin(math.radians(float(_AngSamplingRateDeg)))/4.
      parameters+= \
          ' -perturb ' + str(perturb)

   if (_DoRetricSearchbyTiltAngle):
     parameters+=  \
              ' -min_tilt_angle '      + str(_Tilt0) + \
              ' -max_tilt_angle '      + str(_TiltF)
  
   if (_DoCtfCorrection):
     parameters+=  \
              ' -groups '              + CtfGroupSubsetFileName

   if (_DoParallel):
      parameters = parameters + ' -mpi_job_size ' + str(_MyMpiJobSize)

   launch_job.launch_job('xmipp_angular_project_library',
                         parameters,
                         _mylog,
                         _DoParallel,
                         _MyNumberOfMpiProcesses*_MyNumberOfThreads,
                         1,
                         _MySystemFlavour)


   # Loop over all CTF groups
   # Use reverse order to have same order in add_to docfiles from angular_class_average
   for ii in range(_NumberOfCtfGroups):
      ictf = _NumberOfCtfGroups - ii - 1

      refname          = ProjectLibraryRootName
      if (_DoCtfCorrection):
         CtfGroupName = utils_xmipp.composeFileName(CtfGroupRootName + '_group',ictf+1,'')
         outputname   = ProjMatchRootName + '_' + CtfGroupName 
         CtfGroupName = '../' + CtfGroupDirectory + '/' + CtfGroupName
         inselfile    = CtfGroupName + '.sel'
         inputdocfile = (os.path.basename(inselfile)).replace('.sel','.doc')
         txtfile      = ProjectLibraryRootName + '_sampling.txt'
         if (os.path.exists(txtfile)):
            os.remove(txtfile)
         txtfileb     = utils_xmipp.composeFileName(ProjectLibraryRootName + '_group',ictf+1,'')
         txtfileb     += '_sampling.txt'
         shutil.copy(txtfileb, txtfile)
      else:
         outputname   = ProjMatchRootName
         inputdocfile = _InputDocFileName

      print '*********************************************************************'
      print '* Perform projection matching'
      parameters= ' -i '              + inputdocfile + \
                  ' -o '              + outputname + \
                  ' -ref '            + refname + \
                  ' -Ri '             + str(_Ri)           + \
                  ' -Ro '             + str(_Ro)           + \
                  ' -max_shift '      + str(_MaxChangeOffset) + \
                  ' -search5d_shift ' + str(_Search5DShift) + \
                  ' -search5d_step  ' + str(_Search5DStep) + \
                  ' -mem '            + str(_AvailableMemory) + \
                  ' -sym '            + _SymmetryGroup + 'h'

      if (_DoCtfCorrection and _ReferenceIsCtfCorrected):
         ctffile = CtfGroupName + '.ctf'
         parameters += \
                  ' -pad '            + str(_PaddingFactor) + \
                  ' -ctf '            + ctffile

      if (_DoParallel):
         parameters = parameters + ' -mpi_job_size ' + str(_MyMpiJobSize)

      launch_job.launch_job('xmipp_angular_projection_matching',
                            parameters,
                            _mylog,
                            _DoParallel,
                            _MyNumberOfMpiProcesses*_MyNumberOfThreads,
                            1,
                            _MySystemFlavour)

      # Now make the class averages
      parameters =  ' -i '      + outputname + '.doc'  + \
                    ' -lib '    + ProjectLibraryRootName + '_angles.doc' + \
                    ' -dont_write_selfiles ' + \
                    ' -limit0 ' + str(MinimumCrossCorrelation) + \
                    ' -limitR ' + str(DiscardPercentage)
      if (_DoCtfCorrection):
         # On-the fly apply Wiener-filter correction and add all CTF groups together
         parameters += \
                    ' -wien '             + CtfGroupName + '.wien' + \
                    ' -pad '              + str(_PaddingFactor) + \
                    ' -add_to '           + ProjMatchRootName
      else:
         parameters += \
                    ' -o '                + ProjMatchRootName
      if (_DoAlign2D == '1'):
         parameters += \
                    ' -iter '             + str(_Align2DIterNr) + \
                    ' -Ri '               + str(_Ri)           + \
                    ' -Ro '               + str(_Ro)           + \
                    ' -max_shift '        + str(_MaxChangeOffset) + \
                    ' -max_shift_change ' + str(_Align2dMaxChangeOffset) + \
                    ' -max_psi_change '   + str(_Align2dMaxChangeRot) 
      if (_DoComputeResolution and _DoSplitReferenceImages):
         parameters += \
                    ' -split '

      launch_job.launch_job('xmipp_angular_class_average',
                            parameters,
                            _mylog,
                            _DoParallel,
                            _MyNumberOfMpiProcesses*_MyNumberOfThreads,
                            1,
                            _MySystemFlavour)

      if (_DoAlign2D == '1'):
         outputdocfile =  ProjMatchRootName + '_realigned.doc'
      else:
         outputdocfile =  ProjMatchRootName + '.doc'

      if (_DoCtfCorrection):
         os.remove(outputname + '.doc')
         os.remove(inputdocfile)

   # End loop over all CTF groups

   # Move outputdocfile to standard name
   shutil.move(outputdocfile,_OutputDocFileName)
                
   # Make absolute path so visualization protocol can be run from the same directory
   # Make selfile with reference projections, class averages and realigned averages
   classselfile=selfile.selfile()
   classselfile.read(ProjMatchRootName+'_classes.sel')
   library_sel_file=classselfile.replace_string(ProjMatchRootName+'_class',
                                                        ProjectLibraryRootName);
   before_alignment_sel_file=classselfile.replace_string('.xmp','.ref.xmp');
   before_alignment_sel_file.deactivate_all_images()
   newsel=library_sel_file.intercalate_union_3(before_alignment_sel_file, classselfile)
   compare_sel_file=ProjMatchRootName+'_compare.sel'
   newsel=newsel.make_abspath()
   newsel.write(MultiAlign2dSel)
   # Also make abspath in classes docfile
   newdoc=docfiles.docfile(ProjMatchRootName + '_classes.doc')
   newdoc.make_abspath()
   newdoc.write(ForReconstructionDoc)

   if (_DisplayProjectionMatching):
      command='xmipp_show -sel '+ "../"+'Iter_'+\
                   str(_iteration_number) +'/'+ MultiAlign2dSel +' -w 9 '
      if (_DoAlign2D == '1'):
         command += ' -showall '

      print '*********************************************************************'
      print '* ',command
      _mylog.info(command) 
      os.system(command)

def make_subset_docfiles(_mylog,
                         _InputDocFileName,
                         _NumberOfCtfGroups):

   import os;
   import utils_xmipp
   import launch_job

   # Loop over all CTF groups
   docselfile = []
   for ictf in range(_NumberOfCtfGroups):
      
      CtfGroupName=utils_xmipp.composeFileName(CtfGroupRootName + '_group',ictf+1,'')
      CtfGroupName = '../' + CtfGroupDirectory + '/' + CtfGroupName
      inselfile = CtfGroupName + '.sel'
      inputdocfile = (os.path.basename(inselfile)).replace('.sel','.doc')
      command=' -i   ' + _InputDocFileName + \
              ' -sel ' + inselfile + \
              ' -o   ' + inputdocfile
      print '*********************************************************************'
      launch_job.launch_job("xmipp_docfile_select_subset",
                            command,
                            _mylog,
                            False,1,1,'')
      docselfile.append(inputdocfile+' 1\n')

   # Write the selfile of all these docfiles
   fh = open(CtfGroupSubsetFileName,'w')
   fh.writelines(docselfile)
   fh.close()


#------------------------------------------------------------------------
#execute_reconstruction
#------------------------------------------------------------------------
def execute_reconstruction(_mylog,
                           _ARTReconstructionExtraCommand,
                           _WBPReconstructionExtraCommand,
                           _FourierReconstructionExtraCommand,
                           _iteration_number,
                           _DisplayReconstruction,
                           _DoParallel,
                           _MyNumberOfMpiProcesses,
                           _MyNumberOfThreads,
                           _MySystemFlavour,
                           _MyMpiJobSize,
                           _ReconstructionMethod,
                           _FourierMaxFrequencyOfInterest,
                           _ARTLambda,
                           _SymmetryGroup,
                           _ReconstructedandfilteredVolume,
                           _DoComputeResolution,
                           _DoSplitReferenceImages):

   _mylog.debug("execute_reconstruction")

   import os,shutil,math
   import launch_job

   Outputvolume = _ReconstructedandfilteredVolume

   print '*********************************************************************'
   print '* Reconstruct volume using '
   if _ReconstructionMethod=='wbp':
      Outputvolume = Outputvolume+".vol"
      program = 'xmipp_reconstruct_wbp'
      parameters= ' -i '    + ForReconstructionSel + \
                  ' -doc '  + ForReconstructionDoc + \
                  ' -o '    + Outputvolume + \
                  ' -sym '  + _SymmetryGroup + \
                  ' -weight -use_each_image '
      parameters = parameters + _WBPReconstructionExtraCommand
      _MyNumberOfThreads = 1
              
   elif _ReconstructionMethod=='art':
      program = 'xmipp_reconstruct_art'
      _DoParallel=False
      parameters=' -i '    + ForReconstructionSel + \
                 ' -o '    + Outputvolume + ' ' + \
                 ' -sym '  + _SymmetryGroup + \
		 ' -thr '  + str(_MyNumberOfThreads) + \
                 ' -WLS '
      if len(_ARTLambda)>1:
         parameters = parameters + ' -l '   + _ARTLambda + ' '
      parameters = parameters + _ARTReconstructionExtraCommand
   elif _ReconstructionMethod=='fourier':
      program = 'xmipp_reconstruct_fourier'
      parameters=' -i '    + ForReconstructionSel + \
                 ' -o '    + Outputvolume + '.vol ' + \
                 ' -sym '  + _SymmetryGroup + \
                 ' -thr '  + str(_MyNumberOfThreads) + \
                 ' -weight ' + \
                 ' -max_resolution ' + str(_FourierMaxFrequencyOfInterest)
      if (_DoParallel):
         parameters = parameters + ' -mpi_job_size ' + str(_MyMpiJobSize)
      if (_DoComputeResolution and not _DoSplitReferenceImages):
         parameters = parameters + ' -prepare_fsc ' + Outputvolume + ' '
      parameters = parameters + _FourierReconstructionExtraCommand 
   else:
      _mylog.error("Reconstruction method unknown. Quiting")
      print "Reconstruction method unknown. Quiting"
      exit(1)
    
   launch_job.launch_job(program,
                         parameters,
                         _mylog,
                         _DoParallel,
                         _MyNumberOfMpiProcesses,
                         _MyNumberOfThreads,
                         _MySystemFlavour)

   #_mylog.info(command+ ' ' + parameters)
   if _DisplayReconstruction==True:
      command='xmipp_show -vol '+ Outputvolume + '&'
      print '*********************************************************************'
      print '* ',command
      _mylog.info(command)
      os.system(command)

#------------------------------------------------------------------------
#           execute_resolution(self._SelFileName)
#------------------------------------------------------------------------
def  execute_resolution(_mylog,
                        _ARTReconstructionExtraCommand,
                        _WBPReconstructionExtraCommand,
                        _FourierReconstructionExtraCommand,
                        _ReconstructionMethod,
                        _FourierMaxFrequencyOfInterest,
                        _iteration_number,
                        _DisplayReconstruction,
                        _ResolSam,
                        _DoParallel,
                        _MyNumberOfMpiProcesses,
                        _MyNumberOfThreads,
                        _MySystemFlavour,
                        _MyMpiJobSize,
                        _SymmetryGroup,
                        _DisplayResolution,
                        _ReconstructedVolume,
                        _ARTLambda,
                        _OuterRadius):

    import os,shutil,math

    split_sel_root_name=ProjMatchRootName+'_split'
    Outputvolumes=[]
    Outputvolumes.append(split_sel_root_name+'_1')
    Outputvolumes.append(split_sel_root_name+'_2')
    
    Selfiles=[]
    Selfiles.append(split_sel_root_name+'_1_classes.sel')
    Selfiles.append(split_sel_root_name+'_2_classes.sel')
    Docfiles=[]
    Docfiles.append(split_sel_root_name+'_1_classes.doc')
    Docfiles.append(split_sel_root_name+'_2_classes.doc')
    for i in range(len(Outputvolumes)):
       print '*********************************************************************'
       print '* Reconstruct volume'
       if _ReconstructionMethod=='wbp':
          program = 'xmipp_reconstruct_wbp'
          parameters= ' -i '    + Selfiles[i] + \
                      ' -doc '  + Docfiles[i] + \
                      ' -o '    + Outputvolumes[i] + ".vol" + \
                      ' -sym '  + _SymmetryGroup + \
                      ' -weight -use_each_image '
          parameters = parameters + _WBPReconstructionExtraCommand
          _MyNumberOfThreads = 1
       elif _ReconstructionMethod=='art':
          program = 'xmipp_reconstruct_art'
          _DoParallel=False
          parameters=' -i '    + Selfiles[i] + \
                     ' -o '    + Outputvolumes[i] + \
                     ' -sym '  + _SymmetryGroup + \
		     ' -thr '  + str(_MyNumberOfThreads) + \
                     ' -WLS '
          if len(_ARTLambda)>1:
             parameters = parameters + ' -l '   + _ARTLambda + ' '
          parameters = parameters + _ARTReconstructionExtraCommand
       elif _ReconstructionMethod=='fourier':
          program = 'xmipp_reconstruct_fourier'
          parameters=' -i '    +  Selfiles[i] + \
                     ' -o '    +  Outputvolumes[i] + '.vol ' + \
                     ' -sym '  + _SymmetryGroup + \
		     ' -thr '  + str(_MyNumberOfThreads) + \
                     ' -weight ' + \
                     ' -max_resolution ' + str(_FourierMaxFrequencyOfInterest)
          if (_DoParallel):
             parameters = parameters + ' -mpi_job_size ' + str(_MyMpiJobSize)
          parameters = parameters + _FourierReconstructionExtraCommand
       else:
          _mylog.error("Reconstruction method unknown. Quiting")
          print "Reconstruction method unknown. Quiting"
          exit(1)

       import launch_job
       launch_job.launch_job(program,
                             parameters,
                             _mylog,
                             _DoParallel,
                             _MyNumberOfMpiProcesses,
                             _MyNumberOfThreads,
                             _MySystemFlavour)

    # Prevent high-resolution correlation because of discrete mask from wbp
    innerrad = _OuterRadius - 2
    for i in range(len(Outputvolumes)):
       Outputvolumes[i]+=".vol"
       print '*********************************************************************'
       print '* Applying a soft mask'
       command = " -i " + Outputvolumes[i] + \
                 " -mask  raised_cosine -" + str(innerrad) + \
                 " -" + str(_OuterRadius)
       launch_job.launch_job("xmipp_mask",
                             command,
                             _mylog,
                             False,1,1,_MySystemFlavour)
  
    print '**************************************************************'
    print '* Compute resolution ' 
    command = " -ref " + Outputvolumes[0] +\
              " -i " +Outputvolumes[1]  + ' -sam ' + str(_ResolSam)
    launch_job.launch_job("xmipp_resolution_fsc",
                          command,
                          _mylog,
                          False,1,1,_MySystemFlavour)
    import visualization
    if _DisplayResolution==True:
      plot=visualization.gnuplot()
      plot.plot_xy_file(Outputvolumes[1]+'.frc',
                          Title="Resolution",
                          X_Label="Armstrong^-1",
                          Y_Label="y",
                          X_col=1,
                          Y_col=2)
      print '*********************************************************************'
      print '* plot resolution'
      _mylog.info(" plot resolution")

    # Copy FSC to standard name file
    outputfsc=_ReconstructedVolume.replace(ReconstructedVolume,OutputFsc)
    shutil.copy(Outputvolumes[1]+'.frc',outputfsc) 

    #compute resolution
    resolution_fsc_file = Outputvolumes[1]+'.frc'
    f = open(resolution_fsc_file, 'r')
    #skip first line
    fi=f.readline()
      
    filter_frequence=0. 
    for line in f:
        line = line.strip()
        if not line.startswith('#'):
            mylist = (line.split())
            if( float(mylist[1]) < 0.5):
               break
            else:
              filter_frequence=float(mylist[0])


    f.close()
    print '* maximum resolution (A^-1): ', filter_frequence
    filter_frequence *= _ResolSam
    print '* maximum resolution (px^-1): ', filter_frequence
    return filter_frequence

#------------------------------------------------------------------------
#           filter_at_given_resolution
#------------------------------------------------------------------------
def filter_at_given_resolution(_DoComputeResolution,
                               _mylog, 
                               _iteration_number,
                               _SetResolutiontoZero,
                               _ConstantToAddToFiltration,
                               _filter_frequence,
                               _ReconstructedVolume,
                               _ReconstructedandfilteredVolume,
                               _MySystemFlavour
                               ):

    import os,shutil
    import launch_job
    Inputvolume   =_ReconstructedVolume+'.vol'
    Outputvolume  =_ReconstructedandfilteredVolume+'.vol'
    if (_SetResolutiontoZero):
       filter_in_pixels_at = float(_ConstantToAddToFiltration)
    else:
       filter_in_pixels_at = float(_filter_frequence) +\
                             float(_ConstantToAddToFiltration)
    print '**************************************************************'
    print '* Filter reconstruction ' 
    if (_ConstantToAddToFiltration<0.5 or (not _DoComputeResolution)):
        shutil.copy(Inputvolume,Outputvolume) 
        command ="shutilcopy" + Inputvolume + ' ' + Outputvolume
    else:   
        command = " -i " + Inputvolume +\
                  " -o " + Outputvolume + ' -low_pass ' +\
                  str (filter_in_pixels_at)
        launch_job.launch_job("xmipp_fourier_filter",
                              command,
                              _mylog,
                              False,1,1,_MySystemFlavour)
    return filter_in_pixels_at


#------------------------------------------------------------------------
#create_working directory
#------------------------------------------------------------------------
def execute_cleanup(_mylog,
                    _DeleteClassAverages,
                    _DeleteReferenceProjections):
   import os,glob
   import utils_xmipp
   
   if (_DeleteClassAverages):
      message=' CleanUp: deleting directory '+ ProjMatchDir
      print '* ',message
      _mylog.info(message)
      os.system(' rm -r ' + ProjMatchDir + ' &')

   if (_DeleteReferenceProjections):
      message=' CleanUp: deleting directory '+ LibraryDir
      print '* ',message
      _mylog.info(message)
      os.system(' rm -r ' + LibraryDir + ' &')

def  fill_name_vector(_user_suplied_name,
                      _volume_name_list,
                      _NumberofIterations,
                      _root_name):
     _volume_name_list.append("dummy")
     if (len(_user_suplied_name)>1):
        _volume_name_list.append(_user_suplied_name)
     for _iteration_number in range(1, _NumberofIterations):
         _volume_name_list.append("../"+'Iter_'+\
                                   str(_iteration_number)+'/'+ 'Iter_'+\
                                   str(_iteration_number)+'_'+\
                                   _root_name)
                  
#
# main
#     
if __name__ == '__main__':

    # create rotational_spectra_class object
    # 
   
    #init variables
  my_projmatch=projection_matching_class(
                NumberofIterations,     
                ContinueAtIteration,  
                CleanUpFiles,
                DoMask,   
                DisplayMask,                    
                ReferenceFileName,              
                MaskFileName,                   
                DoProjectionMatching,           
                DisplayProjectionMatching,      
                AngSamplingRateDeg,             
                PerturbProjectionDirections,
                DoRetricSearchbyTiltAngle,      
                Tilt0,                          
                TiltF,                          
                ProjMatchingExtra,              
                MaxChangeOffset,
                MaxChangeInAngles,
                MinimumCrossCorrelation,
                DiscardPercentage,
                DoAlign2D,                      
                InnerRadius,                    
                OuterRadius,                    
                Search5DShift,
                Search5DStep,
                AvailableMemory,
                Align2DIterNr,                  
                Align2dMaxChangeOffset,
                Align2dMaxChangeRot,            
                DisplayReconstruction,
                DisplayResolution,          
                DoReconstruction,
                ReconstructionMethod,
                ARTLambda,
                ARTReconstructionExtraCommand,
                WBPReconstructionExtraCommand,
                FourierReconstructionExtraCommand,
                FourierMaxFrequencyOfInterest,
                DoComputeResolution,
                DoSplitReferenceImages,
                ResolSam,
                SelFileName,                    
                DocFileName,                    
                DoCtfCorrection,
                CTFDatName,
                WienerConstant,
                DoAutoCtfGroup,
                CtfGroupMaxDiff,
                CtfGroupMaxResol,
                SplitDefocusDocFile,
                PaddingFactor,
                DataArePhaseFlipped,
                ReferenceIsCtfCorrected,
                WorkingDir,                  
                ProjectDir,                     
                LogDir,                         
                DoParallel,                     
                NumberOfMpiProcesses,                   
                SystemFlavour,
                MpiJobSize,
                NumberOfThreads,
                SymmetryGroup,                        
                SetResolutiontoZero,
                ConstantToAddToFiltration
                )
