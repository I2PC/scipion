
def expandCommentRun(allowContinue=False):
    list = "Resume, Restart"
    linesStr = '''
#------------------------------------------------------------------------------------------
# {section}{has_question} Comment
#------------------------------------------------------------------------------------------
# Display comment
DisplayComment = False

# {text} Write a comment:
Comment = """Describe your run here..."""
#-----------------------------------------------------------------------------
# {section} Run 
#-----------------------------------------------------------------------------
# RUN name:
""" 
This will identify your protocol run. It need to be unique for each protocol. 
You could have <run1>, <run2> for protocol X, but not two run1 for same protocol. 
This name together with the protocol output folder will determine the working
directory for this run.
"""
RunName = "run_001"

# {list}(%(list)s) Run behavior
""" 
Resume from the last step, restart the whole process 
or continue at a given step or iteration.
"""
Behavior = "Resume"
'''
    if allowContinue:
        list += ", Continue"
        linesStr += '''
# {condition}(Behavior=="Continue") Continue at step:
""" Set to a positive number N to continue the protocol run at step N. """
ContinueAtStep = 1
'''
    return linesStr % locals()

def expandParallel(threads=1, mpi=8, condition="", hours=72, jobsize=0):
    conditionValue = 'True'
    if len(condition) > 0:
        conditionValue = condition
        condition = "{condition}(%s)" % condition
    linesStr = ''' 
#------------------------------------------------------------------------------------------
# {section} %(condition)s Parallelization
#------------------------------------------------------------------------------------------
# {hidden} Parallel Condition
ParallelCondition = "%(conditionValue)s"
'''
    if threads > 0:
        linesStr += '''
# Number of threads
""" 
This option provides shared-memory parallelization on multi-core machines.
It does not require any additional software, other than <Xmipp>
"""
NumberOfThreads = %(threads)d
'''
    if mpi > 0:
        linesStr += '''
# Number of MPI processes
""" 
This option provides the number of independent processes spawned 
in parallel by <mpirun> command in a cluster, usually throught
a queue system. This will require that you have compile <Xmipp>
with <mpi> support.
"""
NumberOfMpi = %(mpi)d
'''
    if jobsize > 0:
        linesStr += '''
        #MPI job size 
"""
Minimum size of jobs in mpi processes. 
Set to 1 for large images (e.g. 500x500)
and to 10 for small images (e.g. 100x100)
"""
MpiJobSize ='%(jobsize)d'
'''
        
    linesStr += '''
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
QueueHours = %(hours)d
''' 
    return linesStr % locals()

def expandExpert():
    return '''    
# {hidden} Show expert options
"""If True, expert options will be displayed """
ShowExpertOptions = False
'''
    
def expandJavaMemory():
    return '''
# {expert} Memory to use (In Gb)
"""
Amount of memory passed to the JVM for the application.
If you have very large micrographs (more than 10k pixels width)
is recommended to use 2Gb or more
"""
Memory = 2
'''
    
def expandParticlesPreprocess(allowFlip):
    linesStr = '''
#-----------------------------------------------------------------------------
# {section} Preprocess
#-----------------------------------------------------------------------------
# Dust removal (Recommended)
""" 
Sets pixels with unusually large values to random values from a Gaussian
with zero-mean and unity-standard deviation. 
"""
DoRemoveDust = True

# {expert}{condition}(DoRemoveDust) Threshold for dust removal:
""" 
Pixels with a signal higher or lower than this value times the standard 
deviation of the image will be affected. For cryo, 3.5 is a good value.
For high-contrast negative stain, the signal itself may be affected so 
that a higher value may be preferable.
"""
DustRemovalThreshold = 3.5
'''
    if allowFlip:
        linesStr += '''
# Phase flipping (Recommended)
""" Use the information from the CTF to compensate for phase reversals."""
DoFlip = True
'''
    linesStr += '''
    
# Invert contrast
""" Invert the contrast if your particles are black over a white background. """
DoInvert = False

# Normalize (Recommended)
""" 
It subtract a ramp in the gray values and normalizes so that in the 
background there is 0 mean and standard deviation 1 """
DoNorm = True

# {expert}{list_combo}(OldXmipp,NewXmipp,Ramp){condition}(DoNorm) Normalization type
"""
OldXmipp (mean(Image)=0, stddev(Image)=1). 
NewXmipp (mean(background)=0, stddev(background)=1), 
Ramp (subtract background+NewXmipp)"""
NormType = "Ramp"

# {expert}{condition}(DoNorm) Background radius
"""
Pixels outside this circle are assumed to be noise and their stddev 
is set to 1. Radius for background circle definition (in pix.).
If this value is 0, then half the box size is used. """
BackGroundRadius = -1
'''
    return linesStr

def expandThreshold():
    return '''
#-----------------------------------------------------------------------------
# {section}{has_question} Threshold 
#-----------------------------------------------------------------------------
# Apply threshold?
DoThreshold = False

# {list_combo}(below, above, abs_below) Select pixels
""" Select those pixels whose value is below, above or whose absolute value is below threshold"""
SelectionMode = 'below' 

# Threshold
""" Threshold value for selecting pixels """
Threshold = 0.0

# {list_combo}(binarize, value, avg) Substitute by
""" Binarize: Selected are set to 0, non-selected to 1; avg: Average of non-selected """
SubstituteBy = 'value'

#{condition}(SubstituteBy=='value') Value
""" Substitute selected pixels by this value """
SubstituteValue = 0.0
'''
    
def expandFilter():
    return '''
#-----------------------------------------------------------------------------
# {section}{has_question} Filter 
#-----------------------------------------------------------------------------
# Apply filters?
DoFilter = False

# Fourier bandpass filter
""" 
You may do a lowpass filter by setting Freq_low to 0. 
You may do a high pass filter by setting Freq_high to 0.5."""
DoFourier = False 

#{condition}(DoFourier){wizard}(wizardChooseBandPassFilter) Freq_low (0<f<0.5)
""" Set it to 0 for low pass filters """
Freq_low = 0.02

#{condition}(DoFourier) {wizard}(wizardChooseBandPassFilter) Freq_high (0<f<0.5)
""" Set it to 0.5 for high pass filters """
Freq_high = 0.35

#{condition}(DoFourier){expert}{wizard}(wizardChooseBandPassFilter) Freq_decay (0<f<0.5)
""" It describes the length of the amplitude decay in a raised cosine """
Freq_decay = 0.02

# Fourier Gaussian
""" Gaussian filter defined in Fourier space"""
DoGaussian = False

#{condition}(DoGaussian){wizard}(wizardChooseGaussianFilter) Frequency sigma
""" Remind that the Fourier frequency is normalized between 0 and 0.5"""
Freq_sigma = 0.04

# Real Gaussian
""" Gaussian filter defined in Real space"""
DoRealGaussian = False

#{condition}(DoRealGaussian){wizard}(wizardChooseRealGaussianFilter) Sigma
""" This sigma is defined in pixel units """
Real_sigma = 2
'''
    
def expandMask():
    return '''
#-----------------------------------------------------------------------------
# {section}{has_question} Mask
#-----------------------------------------------------------------------------
# Apply mask?
""" Apply mask from file """
DoMask = False

# {condition}(DoMask){list_combo}(raised_cosine, circular, binary_file) Mask type
MaskType = "raised_cosine"

# {condition}(DoMask and MaskType!="binary_file"){wizard}(wizardSetBackgroundRadius) Mask radius
MaskRadius = -1

# {condition}(DoMask and MaskType=="raised_cosine") Mask radius outer
MaskRadiusOuter = 2

# {condition}(DoMask and MaskType=="binary_file"){wizard}(wizardDesignMask) Mask file
MaskFile = ""

# {condition}(DoMask){list_combo}(value, min, max, avg) Substitute with
"""Valid values are a value, min, max and avg """
Substitute = "value"

# {condition}(DoMask and Substitute=="value") Substitute value
"""Valid values are a number, min, max and avg """
SubstituteValue = 0
'''
    
def expandResize():
    return '''
#-----------------------------------------------------------------------------
# {section} Resize
#-----------------------------------------------------------------------------
# Resize
""" Change the dimensions of the input images """
DoResize = False

#{condition}(DoResize) New image size
""" New size in pixels of the images """
NewSize = 0

# Crop
"""
This is the desired output size(in pixels) after cropping.
"""
DoCrop = False

# {condition}(DoCrop) Output size:
""" 
In pixels
"""
CropSize = 0
'''

def expandPreprocessFilterMask(allowFlip):
    return expandParticlesPreprocess(allowFlip) + \
        expandFilter() + expandMask()
        
        
from protlib_include_xray import *
from protlib_include_relion import *