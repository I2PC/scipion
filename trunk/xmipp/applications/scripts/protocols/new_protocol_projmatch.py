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
#from XmippData import SingleImgSize
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
SelFileName = 'smallStack.sel'

# {file} {expert} Docfile with the input angles:
""" Do not provide anything if there are no angles yet. 
    In that case, the starting angles will be read from the image headers
    This docfile should be in newXmipp-style format (with filenames as comments)
    Note that all filenames in this docfile should be with absolute paths!
"""
DocFileName = ''

# {file} Initial 3D reference map:
""" Write down the reference/es name. For example "Reference1.vol Reference2.vol"
    specifies two references
"""
ReferenceFileNames = 'volume.spi volum2.spi'

# Working subdirectory: 
""" This directory will be created if it doesn't exist, and will be used to store all output from this run. Don't use the same directory for multiple different runs, instead use a structure like run1, run2 etc. 
"""
WorkingDir = 'ProjMatch/run1'

# Delete working subdirectory if it already exists?
""" Just be careful with this option...
"""
DoDeleteWorkingDir = True

# Number of iterations to perform
NumberofIterations = 4

# Resume at iteration
""" This option may be used to finish a previously performed run.
    Set to 1 to start a new run 
    Note: Do NOT delete working directory if this option is not set to 1
"""
ContinueAtIteration = 1

# {expert} Save disc space by cleaning up intermediate files?
""" Be careful, many options of the visualization protocol will not work anymore, 
    since all class averages, selfiles etc will be deleted.
"""
CleanUpFiles = True

# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project. Often, each data set of a given sample has its own ProjectDir.
"""
ProjectDir = '/home/roberto/tmp/Test'

# {expert} Directory name for logfiles:
LogDir = 'Logs'

#-----------------------------------------------------------------------------
# {section} CTF correction
#-----------------------------------------------------------------------------

# Perform CTF correction?
""" If set to true, a CTF (amplitude and phase) corrected map will be refined,
    and the data will be processed in CTF groups.
    Note that you cannot combine CTF-correction with re-alignment of the classes.
"""
DoCtfCorrection = True

# {file} CTFDat file with CTF data:
""" The input selfile may be a subset of the images in the CTFDat file, but all 
    images in the input selfile must be present in the CTFDat file. This field is 
    obligatory if CTF correction is to be performed. 
    Note that this file should be positioned in the project directory, and that the
    image names and ctf parameter filenames should be in absolute paths.
"""
CTFDatName = 'all_images_new.ctfdat'

# Make CTF groups automatically?
""" Make CTF groups based on a maximum differences at a given resolution limit.
    If this option is set to false, a docfile with the defocus values where to 
    split the images in distinct defocus group has to be provided (see expert option below)
"""
DoAutoCtfGroup = True

# Maximum difference in CTF-values in one group
""" If the difference between the CTF-values up to the resolution limit specified 
    below is larger than the value given here, two images will be placed in 
    distinct CTF groups.
"""
CtfGroupMaxDiff = 0.5

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
SplitDefocusDocFile = ''

# {expert} Padding factor
""" Application of CTFs to reference projections and of Wiener filter to class averages will be done using padded images.
    Use values larger than one to pad the images. Suggestion, use 1 for large image and 2 for small
"""
PaddingFactor = 1.

# {expert} Wiener constant
""" Term that will be added to the denominator of the Wiener filter.
    In theory, this value is the inverse of the signal-to-noise ratio
    If a negative value is taken, the program will use a default value as in FREALIGN 
    (i.e. 10% of average sum terms over entire space) 
    see Grigorieff JSB 157 (2006) pp117-125
"""
WienerConstant = -1

# Images have been phase flipped?
DataArePhaseFlipped = True

# Is the initial reference map CTF (amplitude) corrected?
ReferenceIsCtfCorrected = True

#-----------------------------------------------------------------------------
# {section} Mask
#-----------------------------------------------------------------------------
# Mask reference volume?
""" Masking the reference volume will increase the signal to noise ratio.
    Do not provide a very tight mask.
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mask for details
"""
DoMask = True

# Use a spherical mask?
""" If set to true, provide the radius of the mask in the next input field
    if set to false, provide a binary mask file in the second next input field
"""
DoSphericalMask = True

# Radius of spherical mask
""" This is the radius (in pixels) of the spherical mask 
"""
MaskRadius = 72

# {file} Binary mask file
""" This should be a binary (only 0/1-valued) Xmipp volume of equal dimension as your reference
    The protein region should be white (1) and the solvent should be black (0).
    Note that this entry is only relevant if no spherical mask is used.
"""
MaskFileName = 'mask.vol'

#-----------------------------------------------------------------------------
# {section} Projection Matching
#-----------------------------------------------------------------------------
# Perform projection Matching?
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Projection_matching and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_projection_matching for details
"""
DoProjectionMatching = True

# {expert} Show projection maching library and classes
""" Show average of projections. Do not set this option to true for non-interactive processing (jobs sent to queues)
"""
DisplayProjectionMatching = False

# Inner radius for rotational correlation:
""" In pixels from the image center
"""
InnerRadius = 0

# Outer radius for rotational correlation
""" In pixels from the image center. Use a negative number to use the entire image.
    WARNING: this radius will be use for masking before computing resoution
"""
OuterRadius = 72

# {expert} Available memory to store all references (Gb)
""" This is only for the storage of the references. If your projections do not fit in memory, 
    the projection matching program will run MUCH slower. But, keep in mind that probably 
    some additional memory is needed for the operating system etc.
    Note that the memory per computing node needs to be given. That is, when using threads, 
    this value will be multiplied automatically by the number of (shared-memory) threads.
"""
AvailableMemory = 1

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
AngSamplingRateDeg = '6 4 2 1'

# Angular search range 
"""Maximum change in rot & tilt  (in +/- degrees)
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "1000 1000 10 10 " 
    specifies 4 iterations, the first two set the value to 1000 (no restriction)
    and the last two to 10degrees. An alternative compact notation 
    is ("2x1000 2x10", i.e.,
    2 iterations with value 1000, and 2 with value 10).
    Note: if there are less values than iterations the last value is reused
MaskRadius = 16
    Note: if there are more values than iterations the extra value are ignored
"""
MaxChangeInAngles = '1000 16 12 8 4 2'

# {expert} Perturb projection directions?
""" If set to true, this option will result to a Gaussian perturbation to the 
    evenly sampled projection directions of the reference library. 
    This may serve to decrease the effects of model bias.
"""
PerturbProjectionDirections = False

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
MaxChangeOffset = '1000 10 5'

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
Search5DShift = '4x5 0'

# {expert} Step size for 5D translational search
""" Provide a sequence of numbers (for instance, "2 2 1 1" specifies 4 iterations,
    the first two set the value to 2, then two with 1 pixel.
    An alternative compact notation is ("2x2 2x1", i.e.,
    2 iterations with value 2, and 2 with value 1).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
    
"""
Search5DStep = '2'

# {expert} Restrict tilt angle search?
DoRetricSearchbyTiltAngle = False

# {expert} Lower-value for restricted tilt angle search
Tilt0 = 40

# {expert} Higher-value for restricted tilt angle search
TiltF = 90

# Symmetry group
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry
    for a description of the symmetry groups format
    If no symmetry is present, give c1
"""
SymmetryGroup = 'i3'

# {expert} Symmetry group for Neighbourhood computations
""" If you do not know what this is leave it blank.
    This symmetry will be using for compute neighboring points,
    but not for sampling or reconstruction
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry
    for a description of the symmetry groups format
    If no symmetry is present, give c1
"""
SymmetryGroupNeighbourhood = ''

# {expert} compute only closest neighbor 
""" This option is only relevant if SymmetryGroupNeighbourhood !=''
    If set to True only one neighbor will be computed per sampling point
"""
OnlyWinner = False

# Discard images with ccf below
""" Provide a sequence of numbers (for instance, "0.3 0.3 0.5 0.5" specifies 4 iterations,
    the first two set the value to 0.3, then two with 0.5.
    An alternative compact notation would be ("2x0.3 2x0.5").
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
    Set to -1 to prevent discarding any images
"""
MinimumCrossCorrelation = '-1'

# Discard percentage of images with ccf below
""" Provide a sequence of numbers (for instance, "20 20 10 10" specifies 4 iterations,
    the first two set the value to 20%, then two with 10%
    An alternative compact notation would be ("2x20 2x10").
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
    Set to zero to prevent discarding any images
"""
DiscardPercentage = '10'

# Perform scale search?
""" If true perform scale refinement
"""
DoScale = False

# Step scale factors size
""" Scale step factor size (1 means 0.01 in/de-crements arround 1)
"""
ScaleStep = '1'

# Number of scale steps
""" Number of scale steps.
    With default values (ScaleStep='1' and ScaleNumberOfSteps='3'): 1 +/-0.01 | +/-0.02 | +/-0.03.    
    With values ScaleStep='2' and ScaleNumberOfSteps='4' it performs a scale search over:
     1 +/-0.02 | +/-0.04 | +/-0.06 | +/-0.08.    
"""
ScaleNumberOfSteps = '3'


# {expert} Additional options for Projection_Matching
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Projection_matching and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_projection_matching for details
    try -Ri xx -Ro yy for restricting angular search (xx and yy are
    the particle inner and outter radius)
    
"""
ProjMatchingExtra = ''

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
DoAlign2D = '0'

# {expert} Number of align2d iterations:
""" Use at least 3 iterations
"""
Align2DIterNr = 4

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
Align2dMaxChangeOffset = '2x1000 2x10'

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
Align2dMaxChangeRot = '2x1000 2x20'

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
DoReconstruction = True

# {expert} Display reconstructed volume?
DisplayReconstruction = False

# {list}|fourier|art|wbp| Reconstruction method
""" Choose between wbp, art or fourier
"""
ReconstructionMethod = 'fourier'

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
ARTLambda = '0.2'

# {expert} Additional reconstruction parameters for ART
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Art
        for details
"""
ARTReconstructionExtraCommand = '-k 0.5 -n 10 '

# Initial maximum frequency used by reconstruct fourier
""" This number os only used in the first iteration. 
    From then on, it will be set to resolution computed in the resolution section
"""
FourierMaxFrequencyOfInterest = '0.25'

# {expert} Additional reconstruction parameters for WBP
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Wbp and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_wbp and
        for details
"""
WBPReconstructionExtraCommand = ' '

# {expert} Additional reconstruction parameters for Fourier
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Fourier and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_Fourier and
        for details
    -thr_width 
"""
FourierReconstructionExtraCommand = ' '

#-----------------------------------------------------------------------------
# {section} Compute Resolution
#-----------------------------------------------------------------------------
# Compute resolution?
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Resolution fo details
"""
DoComputeResolution = True

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
"""
DoSplitReferenceImages = True


# Pixel size (in Ang.)
""" This will make that the X-axis in the resolution plots has units 1/Angstrom
"""
ResolSam = 5.6

# {expert} Display resolution?
DisplayResolution = False

#-----------------------------------------------------------------------------
# {section} Postprocessing
#-----------------------------------------------------------------------------
# Low-pass filter the reference?
DoLowPassFilter = True

# Use estimated resolution for low-pass filtering?
"""If set to true, the volume will be filtered at a frecuency equal to
   the  resolution computed with a FSC=0.5 threshold, possibly 
   plus a constant provided by the user in the next input box. 

   If set to false, then the filtration will be made at the constant 
   value provided by the user in the next box (in digital frequency, 
   i.e. pixel-1: minimum 0, maximum 0.5) 
"""
UseFscForFilter = True

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
ConstantToAddToFiltration = '0.1'

# {expert} Center volume
""" Center volume after each iteration """
DoCenterVolume = False

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
DoParallel = False

# Number of MPI processes to use:
NumberOfMpiProcesses = 5

# minumum size of jobs in mpi processe. Set to 1 for large images (e.g. 500x500) and to 10 for small images (e.g. 100x100)
MpiJobSize = '10'

# MPI system Flavour 
""" Depending on your queuing system and your mpi implementation, different mpirun-like commands have to be given.
    Ask the person who installed your xmipp version, which option to use. 
    Or read: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/ParallelPage. The following values are available: 
"""
SystemFlavour = ''

#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript = 'visualize_projmatch.py'

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#-----------------------------------------------------------------------------
#Do not change these variables
ReferenceVolumeName = 'reference_volume.vol'
LibraryDir = "ReferenceLibrary"
ProjectLibraryRootName = LibraryDir + "/ref"
ProjMatchDir = "ProjMatchClasses"
ProjMatchName = 'proj_match'
ProjMatchRootName = ProjMatchDir + "/" + ProjMatchName
ForReconstructionSel = "reconstruction.sel"
ForReconstructionDoc = "reconstruction.doc"
MultiAlign2dSel = "multi_align2d.sel"
DocFileWithOriginalAngles = 'original_angles.doc'
docfile_with_current_angles = 'current_angles.doc'
FilteredReconstruction = "filtered_reconstruction"

ReconstructedVolume = "reconstruction"#
maskReferenceVolume = "mask_reference"#

OutputFsc = "resolution.fsc"
CtfGroupDirectory = "CtfGroups"
CtfGroupRootName = "ctf"
CtfGroupSubsetFileName = "ctf_groups_subset_docfiles.sel"
tableNameRoot = "wrappers"
tableName = tableNameRoot

reconstructedFileNamesIter = []# names for reconstructed volumes
maskedFileNamesIter = []# names masked volumes used as reference
referenceNumber = 1#number of references
createAuxTable = False



import log, logging
import os, sys
#import launch_job
#add to pythonPATH
scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/lib'
sys.path.append(scriptdir) # add default search path
scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/protocols'
sys.path.append(scriptdir)
from pysqlite2 import dbapi2 as sqlite


from xmipp import *
from ProjMatchActionsToBePerformedBeforeLoop import *
from ProjMatchActionsToBePerformedInLoop import *
import pickle

sqlCommand = ""
sqlCommandV = ""

def initDataBase(projectdir, logdir, scriptname, WorkDirectory):
    if logdir[0] == '/':
        LogName = logdir
    else:
        LogName = projectdir + '/' + logdir
    if not LogName[-1] == '/':
        LogName += '/'
    if not os.path.exists(LogName):
        os.makedirs(LogName)
    scriptname = os.path.basename(scriptname)
    LogName += scriptname.replace('.py', '')
    if not (WorkDirectory == "."):
        LogName += '_'
        LogName += os.path.basename(WorkDirectory)
    LogName += '.db'
    conn = sqlite.Connection(LogName)
    # Create table
    #check if table already exists
    #if table exists create an auxiliary table and clean it
    #if user has not requested to start from the beginning
    global tableName
    tableName = tableNameRoot
    sqlCommand = "SELECT count(*) from sqlite_master where tbl_name = ?;"
    cur = conn.cursor()
    cur.execute(sqlCommand, [tableName])
    global createAuxTable
    createAuxTable = cur.fetchone()[0] == 1 and ContinueAtIteration != 1

    if createAuxTable:
        tableName += '_aux'
    sqlCommand = '''create table if not exists ''' + tableName + \
                    '''(id INTEGER PRIMARY KEY,
                     command text, 
                     parameters text,
                     verifyfiles text,
                     init date, 
                     finish date,
                     verified bool)'''

    sqlCommand += ';delete from ' + tableName
    cur.executescript(sqlCommand)
    global sqlCommand
    global sqlCommandV
    sqlCommand = "insert into " + tableName + "(command,parameters)             VALUES (?,?)"
    sqlCommandV = "insert into " + tableName + "(command,parameters,verifyfiles) VALUES (?,?,?)"

    return conn

def compareParameters (conn):
    '''return 0 if new execution of script (tableName) is a subset of and old execution
   this is interesting for continue at iteration option'''
    sqlCommand = '''SELECT count(*) FROM
                          (SELECT * FROM ''' + tableName + ''' 
                           except
                           SELECT * FROM ''' + tableNameRoot + ''' 
                          )'''
    print sqlCommand
    cur = conn.cursor()
    cur.execute(sqlCommand)
    return cur.fetchone()[0]

def actionsToBePerformedBeforeLoopThatDoNotModifyTheFileSystem(_log):
    #1 Convert vectors to list
    from arg import getListFromVector
    global ReferenceFileNames
    ReferenceFileNames = getListFromVector(ReferenceFileNames)
    global referenceNumber
    referenceNumber = len(ReferenceFileNames)

    #name of masked volumes
    global maskedFileNamesIter
    auxList = (referenceNumber + 1) * [None]
    #add dummy name so indexes start a 1
    maskedFileNamesIter.append(auxList)
    for iterN in range(NumberofIterations):
        for refN in range(referenceNumber):
            auxList[refN + 1] = WorkingDir + "/Iter_" + \
                                      str(iterN + 1).zfill(2) + \
                                      '/' + \
                                      maskReferenceVolume + \
                                      "_ref_" + str(refN + 1).zfill(2) + ".vol"
        maskedFileNamesIter.append(list(auxList))

    global reconstructedFileNamesIter
    #add initial reference, useful for mark
    #NOTE THAT INDEXES START AT 1
    reconstructedFileNamesIter.append([""] + ReferenceFileNames)
    for iterN in range(NumberofIterations):
        for refN in range(referenceNumber):
            auxList[refN + 1] = WorkingDir + "/Iter_" + \
                                      str(iterN + 1).zfill(2) + \
                                      '/' + \
                                      ReconstructedVolume + \
                                      "_ref_" + str(refN + 1).zfill(2) + ".vol"
        reconstructedFileNamesIter.append(list(auxList))
    #add initial reference, useful for mark
    reconstructedFileNamesIter.append(ReferenceFileNames)

    #Convert directories/files  to absolute path from projdir
    global CtfGroupDirectory
    CtfGroupDirectory = WorkingDir + '/' + CtfGroupDirectory


def otherActionsToBePerformedBeforeLoop(_log, conn):

    #1Delete working dir
    _Parameters = {
          'ProjectDir':ProjectDir
        , 'WorkingDir':WorkingDir
        , 'DoDeleteWorkingDir':DoDeleteWorkingDir
        , 'ContinueAtIteration':ContinueAtIteration
        }
    parameters = pickle.dumps(_Parameters, 0)
    command = 'deleteWorkingDirectory'
    conn.execute(sqlCommand, [command, parameters])

    #Create directory three
    _Parameters = {
          'ProjectDir':ProjectDir
        , 'WorkingDir':WorkingDir
        , 'NumberofIterations':NumberofIterations
        }
    parameters = pickle.dumps(_Parameters, 0)
    command = 'createRequiredDirectories'
    conn.execute(sqlCommand, [command, parameters])

    #Backup protocol file

    _Parameters = {
          'ProjectDir':ProjectDir
        , 'WorkingDir':WorkingDir
        , 'progName'  :sys.argv[0]
        }
    parameters = pickle.dumps(_Parameters, 0)
    command = 'pm_make_backup_of_script_file'
    conn.execute(sqlCommand, [command, parameters])

    #Check references and projections size match

    _Parameters = {
          'ReferenceFileNames':ReferenceFileNames
        , 'SelFileName':SelFileName
        , 'ContinueAtIteration' :ContinueAtIteration
        }
    parameters = pickle.dumps(_Parameters, 0)
    command = 'checkVolumeProjSize'
    conn.execute(sqlCommand, [command, parameters])

    #Check Option compatibility
    _Parameters = {
          'DoAlign2D':DoAlign2D
        , 'DoCtfCorrection':DoCtfCorrection
        }
    parameters = pickle.dumps(_Parameters, 0)
    command = 'checkOptionsCompatibility'
    conn.execute(sqlCommand, [command, parameters])

    #Init Mask references radius
    _Parameters = {
          'OuterRadius':OuterRadius
        , 'ReferenceFileNames_0':ReferenceFileNames[0]
        , 'SelFileName':SelFileName
        }
    parameters = pickle.dumps(_Parameters, 0)
    command = 'global OuterRadius;OuterRadius = initOuterRadius'
    conn.execute(sqlCommand, [command, parameters])

    #7 make CTF groups
    #too many variables, pass them as a dictionary
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
    parameters = pickle.dumps(_Parameters, 0)
    _VerifyFiles = []
    _VerifyFiles.append(CtfGroupRootName + 'Info.xmd')
    verifyfiles = pickle.dumps(_VerifyFiles, 0)
    command = 'global NumberOfCtfGroups;NumberOfCtfGroups = execute_ctf_groups'
    conn.execute(sqlCommandV, [command, parameters, verifyfiles])

    conn.commit()

def actionsToBePerformedInsideLoop(_log, conn):
    print reconstructedFileNamesIter
    print maskedFileNamesIter
    for iterN in range(ContinueAtIteration, NumberofIterations + 1):
        print "iterN=", iterN
        #############conn.execute(sqlBegin + "MPI_ON" + sqlEnd)
        # Mask reference volume
        for refN in range(1, referenceNumber + 1):
            print "refN=", refN
            # Mask reference volume
            ##############REMOVE SHUTIL.COPY
            _Parameters = {
                                  'DoMask'             : DoMask
                                , 'reconstructedFileName' : reconstructedFileNamesIter[iterN - 1][refN]
                                , 'maskedFileName'     : maskedFileNamesIter[iterN][refN]
                                , 'DoSphericalMask'    : DoSphericalMask
                                , 'maskRadius'         : MaskRadius
                                , 'DoSphericalMask'    : DoSphericalMask
                                , 'userSuppliedMask'   : MaskFileName
                                }
            parameters = pickle.dumps(_Parameters, 0)
            _VerifyFiles = []
            _VerifyFiles.append(maskedFileNamesIter[iterN][refN])
            verifyfiles = pickle.dumps(_VerifyFiles, 0)
            command = "execute_mask"
            conn.execute(sqlCommandV, [command, parameters, verifyfiles])

            #REMOVE
            command = "shutil.copy('%s','%s');dummy" % (ReferenceFileNames[0], reconstructedFileNamesIter[iterN][refN])
            parameters = pickle.dumps("")
            conn.execute(sqlCommand, [command, parameters])

    conn.commit()

def mainLoop(_log, iter):
    conn.row_factory = sqlite.Row
    cur = conn.cursor()
    #check if tableName and tablename_aux are identical if not abort
    if createAuxTable:
        compareParameters(conn)


    cur.execute('''select id, command, parameters 
                   FROM ''' + tableNameRoot +
               ''' WHERE finish IS NOT NULL 
                   ORDER BY id''')

    #Change to Projectdir Execute everything from ProjDir
    os.chdir(ProjectDir)

    for row in cur:
        #print row["command"], row["id"]
        sqlCommand = "update " + tableNameRoot + " set init   = CURRENT_TIMESTAMP where id=%d" % row["id"]
        conn.execute(sqlCommand)
        dict = pickle.loads(str(row["parameters"]))
        print row["command"], _log, dict
        exec (row["command"] + '(_log, dict)')
        sqlCommand = "update " + tableNameRoot + " set finish   = CURRENT_TIMESTAMP where id=%d" % row["id"]
        conn.execute(sqlCommand)
#        >>>>>>>>>>>>>>> verify if adecuate
    conn.commit()
#######
# PROTOCOL STARTS HERE
#######
#create Logging system
# Set up logging
_log = log.init_log_system(ProjectDir,
                            LogDir,
                            sys.argv[0],
                            WorkingDir)

# Uncomment next line to get Debug level logging
_log.setLevel(logging.DEBUG)
_log.debug("Debug level logging enabled")
#init DataBase
conn = initDataBase(ProjectDir,
                    LogDir,
                    sys.argv[0],
                    WorkingDir)

#preprocessing
try:
    actionsToBePerformedBeforeLoopThatDoNotModifyTheFileSystem(_log)
##########################################    otherActionsToBePerformedBeforeLoop(_log, conn)
    actionsToBePerformedInsideLoop(_log, conn)
    mainLoop(_log, iter)
except sqlite.Error, e:
    print "An error occurred:", e.args[0]

                #(in another file) <<<<< define list with actions and parameters
                #                  <<<<< store list in database as it is made
                #                  <<<<< link with restart
#postprocesing
