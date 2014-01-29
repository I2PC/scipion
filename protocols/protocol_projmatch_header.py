#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for projection matching
#
# Example use:
# ./xmipp_protocol_projmatch.py
#
# Authors: Roberto Marabini,
#          Sjors Scheres,    March 2008
#        Rewritten by Roberto Marabini
#

# {begin_of_header}

# {eval} expandCommentRun(allowContinue=True)

#-----------------------------------------------------------------------------
# {section} Input
#-----------------------------------------------------------------------------
# {file}(images*.xmd){validate}(PathExists) Input images:
""" 
Provide a stack or metadata file with the input images.
If you want perform <CTF> correction this file should contains a 
column with <CTFModel> filename. 
The filenames should be relative to the <ProjectDir> 
where you are running the <Protocols>
If you have images outside the <Project> you should import them first.
"""
SelFileName = ''

# {expert} Use initial angles/shifts ? 
""" 
Set to <Yes> if you want to use the initial set of angles/shifts.
This information should be in the images input file. 
(In stack or metadata file).
"""
UseInitialAngles = False

# {file}(*.vol){validate}(PathExists) Initial 3D reference map:
"""
Write down the reference/es name. For example "Reference1.vol Reference2.vol"
specifies two references
"""
ReferenceFileNames = ''

# Number of iterations to perform
NumberOfIterations = 4

# {expert} Clean intermediate files?
""" 
Save disk space by cleaning up intermediate files
mainly volumes
"""
CleanUpFiles = True

#-----------------------------------------------------------------------------
# {section}{has_question} CTF correction
#-----------------------------------------------------------------------------
# Perform CTF correction
""" 
If set to true, a CTF (amplitude and phase) corrected map will be refined,
and the data will be processed in CTF groups.
Note that you cannot combine CTF-correction with re-alignment of the classes.
Remember that CTF information should be provided in the images input file.
"""
DoCtfCorrection = True

# Make CTF groups automatically?
"""
Make CTF groups based on a maximum differences at a given resolution limit.
If this option is set to false, a docfile with the defocus values where to 
split the images in distinct defocus group has to be provided (see expert option below)
"""
DoAutoCtfGroup = True

# Maximum difference in CTF-values in one group
""" If the difference between the CTF-values up to the resolution limit specified 
    below is larger than the value given here, two images will be placed in 
    distinct CTF groups.
"""
CtfGroupMaxDiff = 0.1

# Resolution limit (A) for CTF-differences
""" Maximum resolution where to consider CTF-differences among different groups.
    One should use somewhat higher resolutions than those aimed for in the refinement.
"""
CtfGroupMaxResol = 5.6

# {file} {expert} Docfile with defocus values where to split into groups
""" This field is obligatory if you do not want to make the CTF groups automatically.
    Note that the requested docfile can be made initially with the <xmipp_ctf_group> program,
    and then it can be edited manually to suit your needs. 
"""
SplitDefocusDocFile = ''

# {expert} Padding factor
""" Application of CTFs to reference projections and of Wiener filter to class averages will be done using padded images.
    Use values larger than one to pad the images. Suggestion, use 1 for large image and 2 for small
"""
PaddingFactor = 2

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
    Set to True if reference map has been amplitud corrected
"""
ReferenceIsCtfCorrected =  True

#-----------------------------------------------------------------------------
# {section} {has_question} Mask
#-----------------------------------------------------------------------------
# Mask reference volume
""" Masking the reference volume will increase the signal to noise ratio.
    Do not provide a very tight mask.
    See [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mask] for details
"""
DoMask =True

# Use a spherical mask?
""" If set to true, provide the radius of the mask in the next input field
    if set to false, provide a binary mask file in the second next input field
"""
DoSphericalMask =True

# {condition}(DoSphericalMask){wizard}(wizardSetMaskRadius) Radius of spherical mask
""" This is the radius (in pixels) of the spherical mask 
"""
MaskRadius = -1

# {file} {condition}(not DoSphericalMask)  Binary mask file
""" This should be a binary (only 0/1-valued) Xmipp volume of equal dimension as your reference
    The protein region should be white (1) and the solvent should be black (0).
    Note that this entry is only relevant if no spherical mask is used.
"""
MaskFileName ='mask.vol'

#-----------------------------------------------------------------------------
# {section} Projection Matching
#-----------------------------------------------------------------------------
# {wizard}(wizardSetAlignRadii) Inner radius for rotational correlation:
""" In pixels from the image center
    You may specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "8 8 2 2 " 
    specifies 4 iterations, the first two set the value to 8 
    and the last two to 2. An alternative compact notation 
    is ("2x8 2x0", i.e.,
    2 iterations with value 8, and 2 with value 2).
    <Note>: if there are less values than iterations the last value is reused
    <Note>: if there are more values than iterations the extra value are ignored
"""
InnerRadius = '0'

# {wizard}(wizardSetAlignRadii) Outer radius for rotational correlation
""" In pixels from the image center. Use a negative number to use the entire image.
    <WARNING>: this radius will be use for masking before computing resolution
    You may specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "8 8 2 2 " 
    specifies 4 iterations, the first two set the value to 8 
    and the last two to 2. An alternative compact notation 
    is ("2x8 2x0", i.e.,
    2 iterations with value 8, and 2 with value 2).
    <Note>: if there are less values than iterations the last value is reused
    <Note>: if there are more values than iterations the extra value are ignored
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
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more values than iterations the extra value are ignored
"""
AngSamplingRateDeg='7 5 3 2'

# Angular search range 
"""Maximum change in rot & tilt  (in +/- degrees)
    You may specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "1000 1000 10 10 " 
    specifies 4 iterations, the first two set the value to 1000 (no restriction)
    and the last two to 10degrees. An alternative compact notation 
    is ("2x1000 2x10", i.e.,
    2 iterations with value 1000, and 2 with value 10).
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more values than iterations the extra value are ignored
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
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more values than iterations the extra value are ignored
"""
PerturbProjectionDirections ='0'

# {expert}{list}(fourier, real_space) projection method
"""select projection method, by default real space interpolation is used
"""
ProjectionMethod ='real_space'

# {expert}{condition}(ProjectionMethod=="fourier") Padding factor for projection generation
"""Increase the padding factor will improve projection quality but 
projection generation will be slower. In general padding 1 and spline is OK
"""
PaddingAngularProjection = 1
# {expert} {condition}(ProjectionMethod=="fourier"){list}(nearest,linear,bspline) Interpolation Schema for projection generation
""" Interpolation function.
"""
KernelAngularProjection ='bspline'

# Maximum change in origin offset
""" Maximum allowed change in shift in the 3D+2D searches (in +/- pixels).
    Shifts larger than this value will be reset to (0,0)
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "1000 1000 10 10 " 
    specifies 4 iterations, the first two set the value to 1000 (no restriction)
    and the last two to 10degrees. An alternative compact notation 
    is ("2x1000 2x10", i.e.,
    2 iterations with value 1000, and 2 with value 10).
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more values than iterations the extra value are ignored
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
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more values than iterations the extra value are ignored
    
"""
Search5DShift ='4x5 0'

# {expert} Step size for 5D translational search
""" Provide a sequence of numbers (for instance, "2 2 1 1" specifies 4 iterations,
    the first two set the value to 2, then two with 1 pixel.
    An alternative compact notation is ("2x2 2x1", i.e.,
    2 iterations with value 2, and 2 with value 1).
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more values than iterations the extra value are ignored
    
"""
Search5DStep ='2'

# {expert} Restrict tilt angle search?
DoRestricSearchbyTiltAngle = False

# {expert} {condition}(DoRestricSearchbyTiltAngle) Lower-value for restricted tilt angle search
Tilt0 = -91

# {expert} {condition}(DoRestricSearchbyTiltAngle) Higher-value for restricted tilt angle search
TiltF = 91

# Symmetry group
""" See [http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry]
    for a description of the symmetry groups format
    If no symmetry is present, give c1
"""
SymmetryGroup ='c1'

# {expert} Symmetry group for Neighbourhood computations
""" If you do not know what this is leave it blank.
    This symmetry will be using for compute neighboring points,
    but not for sampling or reconstruction
    See [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry]
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
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more values than iterations the extra value are ignored
"""
OnlyWinner ='0'

# {list}(none, maxCC, percentage, classPercentage) Discard images?
""" Choose between none, maxCC, percentage, classPercentage. 
    none : No images will be discarded.
    maxCC  : Minimum Cross Correlation, discard images with CC below a fixed value.
    percentage : Discard percentage of images with less CC.
    classPercentage: Discard percentage of images in each projection direction with less CC.
    Value of each option is set below.
"""
DiscardImages = 'none'

# {condition}(DiscardImages=="maxCC") with CC below
""" 
    Discard images with cross-correlation (CC) below this value.
    Provide a sequence of numbers (for instance, "0.3 0.3 0.5 0.5" specifies 4 iterations,
    the first two set the value to 0.3, then two with 0.5.
    An alternative compact notation would be ("2x0.3 2x0.5").
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more values than iterations the extra value are ignored
"""
MinimumCrossCorrelation = '0.1'

# {condition}(DiscardImages=="percentage") images percent with less CC
""" 
    Discard this percentage of images with less cross-correlation (CC)
    Provide a sequence of numbers (for instance, "20 20 10 10" specifies 4 iterations,
    the first two set the value to 20%, then two with 10%
    An alternative compact notation would be ("2x20 2x10").
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more values than iterations the extra value are ignored
    Set to zero to prevent discarding any images
"""
DiscardPercentage ='10'

# {condition}(DiscardImages=="classPercentage") images percent in each class with less CC
""" 
    Discard this percentage of images in each class(projection direction)
    with less cross-correlation (CC)    
    Provide a sequence of numbers (for instance, "20 20 10 10" specifies 4 iterations,
    the first two set the value to 20%, then two with 10%
    An alternative compact notation would be ("2x20 2x10").
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more values than iterations the extra value are ignored
    Set to zero to prevent discarding any images
"""
DiscardPercentagePerClass ='10'


# Perform scale search?
""" 
    If true perform scale refinement. (Under development)
"""
DoScale = False

# {condition}(DoScale) Step scale factors size
""" Scale step factor size (1 means 0.01 in/de-crements arround 1). MAybe different for the different iterations
"""
ScaleStep ='1'

# {condition}(DoScale) Number of scale steps
""" 
    Number of scale steps.
    With default values (ScaleStep='1' and ScaleNumberOfSteps='3'): 1 +/-0.01 | +/-0.02 | +/-0.03.    
    With values ScaleStep='2' and ScaleNumberOfSteps='4' it performs a scale search over:
     1 +/-0.02 | +/-0.04 | +/-0.06 | +/-0.08.    
    In general scale correction should only be applied to the last iteration. Do not use it unless
    your data is fairly well aligned.

"""
ScaleNumberOfSteps ='3'


# {expert} Additional options for Projection_Matching
""" For details see:
    [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Projection_matching] and
    [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_projection_matching]
    try -Ri xx -Ro yy for restricting angular search (xx and yy are
    the particle inner and outter radius)
    
"""
ProjMatchingExtra =''

# {expert}Save images assigned to each class?
""" If true, save images assigned to each class to a metadata file
    Be aware that for a very fine angular sampling it can be time consuming.
    Not a very parameter any longer since you can get this information in 
    visualize results
"""
DoSaveImagesAssignedToClasses = False

#-----------------------------------------------------------------------------
# {section}{expert}{has_question} 2D re-alignment of classes
#-----------------------------------------------------------------------------
# Perform 2D re-alignment
PerformAlign2D = False

# Perform 2D re-alignment of classes?
""" After performing a 3D projection matching iteration, each of the
    subsets of images assigned to one of the library projections is
    re-aligned using a 2D-alignment protocol.
    This may serve to remove model bias.
    For details see:
    [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Align2d]
    Note that you cannot combine this option with CTF-correction!
    You may specify this option for each iteration. 
    This can be done by a sequence of 0 or 1 numbers (for instance, "1 1 0 0" 
    specifies 4 iterations, the first two applied alig2d while the last 2
    dont. an alternative compact notation is 
    is ("2x1 2x0", i.e.,
    2 iterations with value 1, and 2 with value 0).
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more values than iterations the extra value are ignored
    <IMPORTANT:> if you set this variable to 0 the output  of the projection
    muching step will be copied as output of align2d
"""
DoAlign2D ='0'

# Number of align2d iterations:
""" Use at least 3 iterations
    The number of align iteration may change in each projection matching iteration
    Ffor instance, "4 4 3 3 " 
    specifies 4 alig2d iterations in the first projection matching iteration 
    and  two 3 alig2d iteration in the last 2 projection matching iterations.
     An alternative compact notation 
    is ("2x4 2x3", i.e.,
    2 iterations with value 4, and 2 with value 3).
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more values than iterations the extra value are ignored
"""
Align2DIterNr ='4'

# Maximum change in origin offset (+/- pixels)
"""Maximum change in shift  (+/- pixels)
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "1000 1000 10 10 " 
    specifies 4 iterations, the first two set the value to 1000 (no restriction)
    and the last two to 10degrees. An alternative compact notation 
    is ("2x1000 2x10", i.e.,
    2 iterations with value 1000, and 2 with value 10).
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more values than iterations the extra value are ignored
"""
Align2dMaxChangeOffset ='2x1000 2x10'

# Maximum change in rotation (+/- degrees)
"""Maximum change in shift  (+/- pixels)
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "1000 1000 10 10 " 
    specifies 4 iterations, the first two set the value to 1000 (no restriction)
    and the last two to 10degrees. An alternative compact notation 
    is ("2x1000 2x10", i.e.,
    2 iterations with value 1000, and 2 with value 10).
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more values than iterations the extra value are ignored
"""
Align2dMaxChangeRot ='2x1000 2x20'

#-----------------------------------------------------------------------------
# {section} 3D Reconstruction
#-----------------------------------------------------------------------------

# {list}(fourier, art, wbp) Reconstruction method
""" Choose between wbp, art or fourier
"""
ReconstructionMethod ='fourier'

# {expert}{condition}(ReconstructionMethod=="art") Values of lambda for art
""" <IMPORTANT>: ou must specify a value of lambda for each iteration even
    if art has not been selected.
    <IMPORTANT:> NOte that we are using the WLS version of ART that 
    uses geater lambdas than the plain art.
    See for details:
    [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Art]
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, ".1 .1 .3 .3" 
    specifies 4 iterations, the first two set the value to 0.1 
    (no restriction)
    and the last  two to .3. An alternative compact notation 
    is ("2x.1 2x.3").
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more values than iterations the extra value are ignored
"""
ARTLambda ='0.2'

# {expert}{condition}(ReconstructionMethod=="art") Additional reconstruction parameters for ART
""" For details see:
    [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Art]
"""
ARTReconstructionExtraCommand ='-k 0.5 -n 10 '

# {condition}(ReconstructionMethod=="fourier") Initial maximum frequency
""" This number is only used in the first iteration. 
    From then on, it will be set to resolution computed in the resolution section
"""
FourierMaxFrequencyOfInterest =0.25

# {expert}{condition}(ReconstructionMethod=="wbp") Additional reconstruction parameters for WBP
""" For details see:
    [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Wbp] and
    [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_wbp]
"""
WBPReconstructionExtraCommand =''

# {expert} {condition}(ReconstructionMethod=="fourier")Additional reconstruction parameters for Fourier
""" For details see:
    [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Fourier] and
    [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_Fourier] and
    -thr_width 
"""
FourierReconstructionExtraCommand =''

#-----------------------------------------------------------------------------
# {section} Compute Resolution
#-----------------------------------------------------------------------------
# Compute resolution?
""" For details see:
    [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Resolution]
    Set to 1 to compute resolution and to 0 if you do not want to compute it.

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
   <IMPORTANT:> the second option has ONLY been implemented for FOURIER
   reconstruction method. Other reconstruction methods require this
   flag to be set to True
    You may specify this option for each iteration. 
    This can be done by a sequence of 0 or 1 numbers (for instance, "1 1 0 0" 
    specifies 4 iterations, the first two split the images   while the last 2
    don't. an alternative compact notation is 
    is ("2x1 2x0", i.e.,
    2 iterations with value 1, and 2 with value 0).
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more vapplications/scripts/protocols/new_protocol_projmatch.pyalues than iterations the extra value are ignored
"""
DoSplitReferenceImages ="1"


## sampling is now in acquisition info file
## Pixel size (in Ang.)
#""" This will make that the X-axis in the resolution plots has units 1/Angstrom
#"""
#ResolSam=5.6

#-----------------------------------------------------------------------------
# {section} Postprocessing
#-----------------------------------------------------------------------------
# Low-pass filter the reference?
DoLowPassFilter = True

# {condition}(DoLowPassFilter) Use estimated resolution for low-pass filtering?
"""If set to true, the volume will be filtered at a frecuency equal to
   the  resolution computed with a FSC=0.5 threshold, possibly 
   plus a constant provided by the user in the next input box. 

   If set to false, then the filtration will be made at the constant 
   value provided by the user in the next box (in digital frequency, 
   i.e. pixel-1: minimum 0, maximum 0.5) 
"""
UseFscForFilter = True

# {condition}(DoLowPassFilter) Constant to be added to the estimated resolution
""" The meaning of this field depends on the previous flag.
    If set to true, then the volume will be filtered at a frecuency equal to
    the  resolution computed with resolution_fsc (FSC=0.5) plus the value 
    provided in this field 
    If set to false, the volume will be filtered at the resolution
    provided in this field 
    This value is in digital frequency, or pixel^-1: minimum 0, maximum 0.5

    If you detect correlation between noisy regions decrease this value 
    (even to negative values)
    
    You can specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, ".15 .15 .1 .1" 
    specifies 4 iterations, the first two set the constant to .15
    and the last two to 0.1. An alternative compact notation 
    is ("2x.15 2x0.1", i.e.,
    4 iterations with value 0.15, and three with value .1).
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more values than iterations the extra value are ignored
"""
ConstantToAddToFiltration ='0.05'

# Constant to be added to the reconstruction maximum frequency
""" The meaning of this field depends on the UseFscForFilter flag.
    If set to true, then the volume will be reconstructed up to the frequency equal to
    the resolution computed with resolution_fsc (FSC=0.5) plus the value 
    provided in this field 
    If set to false, the volume will be reconstructed up to the resolution
    provided in this field 
    This value is in digital frequency, or pixel^-1: minimum 0, maximum 0.5

    You can specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, ".15 .15 .1 .1" 
    specifies 4 iterations, the first two set the constant to .15
    and the last two to 0.1. An alternative compact notation 
    is ("2x.15 2x0.1", i.e.,
    4 iterations with value 0.15, and three with value .1).
    <Note:> if there are less values than iterations the last value is reused
    <Note:> if there are more values than iterations the extra value are ignored
"""
ConstantToAddToMaxReconstructionFrequency ='0.1'

# {eval} expandParallel(jobsize=1)

#------------------------------------------------------------------------------------------------
# {section}{visualize} Visualization
#------------------------------------------------------------------------------------------------

# Show results for iterations
""" You may specify more than one iteration here 
    This can be done by a sequence of numbers (for instance, "2 8" 
    specifies iteration 2 and 8 (but not 3, 4, 5, 6 and 7)
"""
DisplayIterationsNo='1 2 3'

# Show results for reference 3D volumes
""" 
   If you want two see the reference volume 2 and 5 write
   2 5
"""
DisplayRef3DNo='1'

# {expert} Width of projection galleries
""" 
    Usually a multiple of 2 is the right value. -1 => authomatic
"""
MatrixWidth=-1

#------------------------------------------------------------------------------------------------
# {section}{visualize} Volumes display
#------------------------------------------------------------------------------------------------

# {list}(x, y, z, surface) Display volumes as slices or surface rendering
""" x -> Visualize volumes in slices along x
    y -> Visualize volumes in slices along y
    z -> Visualize volumes in slices along z
    For surface rendering to work, you need to have chimera installed!
"""
DisplayVolumeSlicesAlong='z'

# {view} Display reference volume
""" Volume after filtration and masking
"""
DisplayReference = False
# {view} Display reconstructed volume
""" Volume as given by the reconstruction algorithm
"""
DisplayReconstruction=False
# {view} Display reconstructed volume after filtration
""" Volume after filtration
"""
DisplayFilteredReconstruction=True

# {view} Display a b_factor corrected volume
""" This utility boost up the high frequencies. Do not use the automated 
    mode [default] for maps with resolutions lower than 12-15 Angstroms.
    It does not make sense to apply the Bfactor to the firsts iterations
    see http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Correct_bfactor.
    NOTE: bfactor will be applied ONLY to the reconstructed volumes NOT to
    the reference ones
"""
DisplayBFactorCorrectedVolume=False

#### this information is in acquisition info
#### {condition}(DisplayBFactorCorrectedVolume) Sampling rate
#### SamplingRate=5.6

#{condition}(DisplayBFactorCorrectedVolume) Maximum resolution to apply B-factor (in Angstrom)
MaxRes=12

# {expert} User defined flags for the correct_bfactor program 
""" See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Correct_bfactor
    for details. DEFAULT behaviour is --auto
"""
CorrectBfactorExtraCommand='--auto'

#------------------------------------------------------------------------------------------------
# {section}{visualize} Library, classes and images
#------------------------------------------------------------------------------------------------

# {view} Display library?
DisplayProjectionMatchingLibrary=False

# {view} Display classes?
DisplayProjectionMatchingClasses=False

# {view} Display library and classes in a single image?
DisplayProjectionMatchingLibraryAndClasses=False

# {view} Display library and experimental images in a single image?
DisplayProjectionMatchingLibraryAndImages=False

# {view} Display discarded images?
DisplayDiscardedImages=False

#-----------------------------------------------------------------------------
# {section}{visualize} Convergence 
#-----------------------------------------------------------------------------
# {view} Plot histogram with angular/shift changes
""" Plot histogram with angular changes from one iteration to next. 
Iteration 0 -> initial values
"""
PlotHistogramAngularMovement=False

# {expert}{condition}(PlotHistogramAngularMovement) Number of bins (for histogram)
""" Number of bins in histograms
"""
NumberOfBins=50

# {condition}(PlotHistogramAngularMovement) Use Psi to compute angular distances
""" Use Psi
"""
UsePsi = False

# {condition}(PlotHistogramAngularMovement) Sort angular metadata and save it.
""" Sort metadata with  experimental images a save it. 
Sorting made by angle
"""
AngleSort = False 
# {condition}(PlotHistogramAngularMovement) Sort shift metadata and save it.
""" Sort metadata with  experimental images a save it. 
Sorting made by shift
"""
ShiftSort = False 

# {condition}(PlotHistogramAngularMovement) save worst particles (percentage) 
"""save a metadata file with the particles 
than move more than a 100-xx percent in all the iteration selected.
If negative file will not be saved.
"""
Percentage = -1 

#------------------------------------------------------------------------------------------------
# {section}{visualize} Angular distribution and resolution plots
#------------------------------------------------------------------------------------------------

# {view} Display angular distribution?
DisplayAngularDistribution=True

# {condition}(DisplayAngularDistribution) {list} (2D, 3D) Display Angular distribution with
""" 2D option uses matplotlib while 3D chimera
"""
DisplayAngularDistributionWith='2D'

# {view} Display resolution plots (FSC)
DisplayResolutionPlots=True

# {expert} Display a threshold in resolution plots (FSC)
ResolutionThreshold=0.5


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#-----------------------------------------------------------------------------
       
from protocol_projmatch import *

if __name__ == '__main__':
    protocolMain(ProtProjMatch)
