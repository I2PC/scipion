# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-package implement projection matching using xmipp 3.1
"""
"""TODO:

1) expert for sections
2) get set list parameters
3) Data type for symmetry with verify
"""
#/home/josem/work/development/testXmipp/tmpProject/BPV_Project
#execute: pw_project.py TestXmippWorkflow

from pyworkflow.em import *  
from pyworkflow.utils import *  
import xmipp, xmipp3

        
class XmippProtProjMatch(xmipp3.XmippProtocol, ProtRefine3D, ProtClassify3D):
    """ Protocol for Xmipp-based ProjMatch/MLF3D classification. """
    _label = 'projection matching'
    
    # Reconstruction method constants
    FOURIER = 0
    WLSART = 1

    def _defineParams(self, form):
        form.addSection(label='Input')
        
        #SelFileName
        form.addParam('inputParticles', PointerParam, label="Input particles", important=True, 
                      pointerClass='SetOfParticles', 
                      help='Select the input particles. \n '
                           'If you want perform *CTF* correction the input particles \n '
                           'should have information about the CTF (hasCTF=True)')  
        form.addParam('useInitialAngles', BooleanParam, default=False,
                      label="Use initial angles/shifts ? ", 
                      help='Set to *Yes* if you want to use the projection assignment (angles/shifts) \n '
                      'associated with the input particles (hasProjectionAssigment=True)')
        # ReferenceFileNames      
        form.addParam('input3DReferences', PointerParam,
                      pointerClass='SetOfVolumes',
                      label='Initial 3D reference volumes', 
                      help='Initial 3D density maps with the same dimensions as your particles. \n '
                           'For example: reference1.vol reference2.vol \n '
                           'specifies two references.')
        form.addParam('numberOfIterations', IntParam, default=4,
                      label='Number of iterations',
                      help='Number of iterations to perform.')
        form.addParam('cleanUpFiles', BooleanParam, default=False,
                      label="Clean up intermediate files?",  expertLevel=LEVEL_EXPERT,
                      help='Save disc space by cleaning up intermediate files. \n '
                           'Be careful, many options of the visualization protocol will not work anymore, \n '
                           'since all class averages, selfiles etc will be deleted. ')
        
        form.addSection(label='CTF correction')
        
        form.addParam('doCTFCorrection', BooleanParam, default=True,
                      label="Perform CTF correction?", 
                      help='If set to true, a CTF (amplitude and phase) corrected map will be refined, \n '
                           'and the data will be processed in CTF groups. \n '
                           'Note that you cannot combine CTF-correction with re-alignment of the classes. \n '
                           'Remember that CTF information should be provided in the images input file. \n ')    
              
        form.addParam('doAutoCTFGroup', BooleanParam, default=True, condition='doCTFCorrection',
                      label="Make CTF groups automatically?", 
                      help='Make CTF groups based on a maximum differences at a given resolution limit. \n '
                           'If this option is set to false, a docfile with the defocus values where to  \n '
                           'split the images in distinct defocus group has to be provided (see expert option below) \n ')             
              
        form.addParam('ctfGroupMaxDiff', FloatParam, default=0.1, condition='doCTFCorrection and doAutoCTFGroup',
                      label='Maximum difference for grouping', validators=[Positive],
                      help='If the difference between the CTF-values up to the resolution limit specified \n '
                      'below is larger than the value given here, two images will be placed in \n '
                      'distinct CTF groups.')          
        
        form.addParam('ctfGroupMaxResol', FloatParam, default=5.6, condition='doCTFCorrection and doAutoCTFGroup',
                      label='Resolution limit (Ang) for grouping', validators=[Positive],
                      help='Maximum resolution where to consider CTF-differences among different groups. \n '
                            'One should use somewhat higher resolutions than those aimed for in the refinement.')       
        
        # SplitDefocusDocFile
        form.addParam('setOfDefocus', StringParam,
                      label='Set of defocus', default='', condition='doCTFCorrection and not doAutoCTFGroup',
                      help='Set with defocus values where to split into groups. \n '
                           'This field is compulsory if you do not want to make the CTF groups automatically. \n '
                           'Note that the requested docfile can be made initially with the *xmipp_ctf_group* program, \n '
                           'and then it can be edited manually to suit your needs.')
        
        form.addParam('paddingFactor', FloatParam, default=2, condition='doCTFCorrection',
                      label='Padding factor', validators=[GE(1)],
                      help='Application of CTFs to reference projections and of Wiener filter \n '
                            'to class averages will be done using padded images. \n '
                            'Use values larger than one to pad the images.')        
        
        form.addParam('wienerConstant', FloatParam, default=-1, condition='doCTFCorrection',
                      label='Wiener constant',  expertLevel=LEVEL_EXPERT,
                      help='Term that will be added to the denominator of the Wiener filter. \n '
                            'In theory, this value is the inverse of the signal-to-noise ratio \n '
                            'If a negative value is taken, the program will use a default value as in FREALIGN \n '
                            '(i.e. 10% of average sum terms over entire space)  \n '
                            'see Grigorieff JSB 157 (2006) pp117-125')   
        
        #TODO: Use common mask parameters
        
        form.addSection(label='Mask')
        
        # doMask, doSphericalMask now merged into maskType
        
        form.addParam('maskType', EnumParam, choices=['None', 'circular', 'binary file'], default=xmipp3.MASK_CIRCULAR, 
                      label="Mask reference volumes", display=EnumParam.DISPLAY_COMBO,
                      help='Masking the reference volume will increase the signal to noise ratio. \n '
                           'Do not provide a very tight mask. \n ')
        
        form.addParam('maskRadius', IntParam, default=-1, condition='maskType == 1',
                      label='Radius of spherical mask (pix)',
                      help='This is the radius (in pixels) of the spherical mask ')       
        
        form.addParam('maskFile', StringParam, default='maks.vol', 
                      label='Binary mask file', condition='maskType == 2',
                      help='The mask file should have the same dimensions as your input particles. \n '
                           'The protein region should be 1 and the solvent should be 0.')  
        
        # DataArePhaseFlipped , now taken from inputParticles.isPhaseFlipped()
        # ReferenceIsCtfCorrected, now taken from input3DReferences.isAmplitudeCorrected()
        
        form.addSection(label='Projection Matching')
        
        form.addParam('innerRadius', NumericListParam, default='0', 
                      label='Inner radius for rotational correlation:', 
                      help=""" In pixels from the image center
    You may specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "8 8 2 2 " 
    specifies 4 iterations, the first two set the value to 8 
    and the last two to 2. An alternative compact notation 
    is ("2x8 2x0", i.e.,
    2 iterations with value 8, and 2 with value 2).
    *Note:* if there are less values than iterations the last value is reused
    *Note:* if there are more values than iterations the extra value are ignored
""")
              
        form.addParam('outerRadius', NumericListParam, default='64', 
                      label='Outer radius for rotational correlation', 
                      help=""" In pixels from the image center. Use a negative number to use the entire image.
*WARNING*: this radius will be use for masking before computing resolution
You may specify this option for each iteration. 
This can be done by a sequence of numbers (for instance, "8 8 2 2 " 
specifies 4 iterations, the first two set the value to 8 
and the last two to 2. An alternative compact notation 
is ("2x8 2x0", i.e.,
2 iterations with value 8, and 2 with value 2).
*Note:* if there are less values than iterations the last value is reused
*Note:* if there are more values than iterations the extra value are ignored
""")        
        
        form.addParam('availableMemory', IntParam, default=2, expertLevel=LEVEL_ADVANCED, 
                      label='Available memory to store all references (Gb)',
                      help=""" This is only for the storage of the references. If your projections do not fit in memory, 
the projection matching program will run MUCH slower. But, keep in mind that probably 
some additional memory is needed for the operating system etc.
Note that the memory per computing node needs to be given. That is, when using threads, 
this value will be multiplied automatically by the number of (shared-memory) threads.
""")
        
        form.addParam('angSamplingRateDeg', NumericListParam, default='7 5 3 2', 
                      label='Angular sampling rate (deg)',
                      help=""" Angular distance (in degrees) between neighboring projection  points
You may specify this option for each iteration. 
This can be done by a sequence of numbers (for instance, "8 8 2 2 " 
specifies 4 iterations, the first two set the value to 8 
and the last two to 2. An alternative compact notation 
is ("2x8 2x0", i.e.,
2 iterations with value 8, and 2 with value 2).
*Note:* if there are less values than iterations the last value is reused
*Note:* if there are more values than iterations the extra value are ignored
""")
        
#        form.addParam('maxChangeInAngles', NumericListParam, default='1000 10 4 2', 
#                      label='Angular search range (deg)',
#                      help=""" Maximum change in rot & tilt  (in +/- degrees)
#    You may specify this option for each iteration. 
#    This can be done by a sequence of numbers (for instance, "1000 1000 10 10 " 
#    specifies 4 iterations, the first two set the value to 1000 (no restriction)
#    and the last two to 10degrees. An alternative compact notation 
#    is ("2x1000 2x10", i.e.,
#    2 iterations with value 1000, and 2 with value 10).
#    <Note:> if there are less values than iterations the last value is reused
#    <Note:> if there are more values than iterations the extra value are ignored
#""")        
        
        form.addParam('perturbProjectionDirections', NumericListParam, default='0', 
                      label='Perturb projection directions?', expertLevel=LEVEL_EXPERT,
                      help=""" If set to 1, this option will result to a Gaussian perturbation to the 
evenly sampled projection directions of the reference library. 
This may serve to decrease the effects of model bias.
You may specify this option for each iteration. 
This can be done by a sequence of numbers (for instance, "1 1 0" 
specifies 3 iterations, the first two set the value to 1 
and the last to 0. An alternative compact notation 
is ("2x1 0", i.e.,
2 iterations with value 1, and 1 with value 0).
*Note:* if there are less values than iterations the last value is reused
*Note:* if there are more values than iterations the extra value are ignored
""")   
        
        # Changed from String to Int 
        form.addParam('projectionMethod', EnumParam, choices=['fourier', 'real space'], 
                      default=xmipp3.PROJECT_FOURIER, expertLevel=LEVEL_EXPERT, 
                      label="Projection method", display=EnumParam.DISPLAY_COMBO,
                      help='select projection method, by default Fourier with padding 1 and interpolation bspline')        

        
        form.addParam('paddingAngularProjection', FloatParam, default=1, expertLevel=LEVEL_EXPERT,  
                      condition='projectionMethod == %d' % xmipp3.PROJECT_FOURIER,
                      label='Padding factor for projection', validators=[GE(1)],
                      help="""Increase the padding factor will improve projection quality but 
projection generation will be slower. In general padding 1 and spline is OK
""")       
        # Changed from String to Int 
        form.addParam('kernelAngularProjection', EnumParam, choices=['neareast', 'linear', 'bspline'],
                      default=xmipp3.KERNEL_BSPLINE, expertLevel=LEVEL_EXPERT,  
                      condition='projectionMethod == %d' % xmipp3.PROJECT_FOURIER,
                      label='Interpolation kernel for projection', 
                      help=""" Interpolation kernel for the generation of projections.
""")
        
        form.addParam('maxChangeOffset', NumericListParam, default='1000 10 5', 
                      label='Perturb projection directions?', expertLevel=LEVEL_EXPERT,
                      help=""" If set to 1, this option will result to a Gaussian perturbation to the 
evenly sampled projection directions of the reference library. 
This may serve to decrease the effects of model bias.
You may specify this option for each iteration. 
This can be done by a sequence of numbers (for instance, "1 1 0" 
specifies 3 iterations, the first two set the value to 1 
and the last to 0. An alternative compact notation 
is ("2x1 0", i.e.,
2 iterations with value 1, and 1 with value 0).
*Note:* if there are less values than iterations the last value is reused
*Note:* if there are more values than iterations the extra value are ignored
""")          

        form.addParam('search5DShift', NumericListParam, default='4x5 0', 
                      label='Search range for 5D translational search',
                      help=""" Give search range from the image center for 5D searches (in +/- pixels).
Values larger than 0 will results in 5D searches (which may be CPU-intensive)
Give 0 for conventional 3D+2D searches. 
Note that after the 5D search, for the optimal angles always 
a 2D exhaustive search is performed anyway (making it ~5D+2D)
Provide a sequence of numbers (for instance, "5 5 3 0" specifies 4 iterations,
the first two set the value to 5, then one with 3, resp 0 pixels.
An alternative compact notation is ("3x5 2x3 0", i.e.,
3 iterations with value 5, and 2 with value 3 and the rest with 0).
*Note:* if there are less values than iterations the last value is reused
*Note:* if there are more values than iterations the extra value are ignored

""")  
        form.addParam('search5DStep', NumericListParam, default='2', 
                      label='Step size for 5D translational search', expertLevel=LEVEL_EXPERT,
                      help="""" Provide a sequence of numbers (for instance, "2 2 1 1" specifies 4 iterations,
    the first two set the value to 2, then two with 1 pixel.
    An alternative compact notation is ("2x2 2x1", i.e.,
    2 iterations with value 2, and 2 with value 1).
    *Note:* if there are less values than iterations the last value is reused
    *Note:* if there are more values than iterations the extra value are ignored
""")          

        form.addParam('doRestricSearchbyTiltAngle', BooleanParam, default=False, expertLevel=LEVEL_EXPERT,
                      label="Restrict tilt angle search?", 
                      help ='Restrict tilt angle search \n ')             

        form.addParam('tilt0', FloatParam, default=-91., condition='doRestricSearchbyTiltAngle',
                      label="Lower-value for restricted tilt angle search", 
                      help ='Lower-value for restricted tilt angle search \n ')             

        form.addParam('tiltF', FloatParam, default=-91., condition='doRestricSearchbyTiltAngle',
                      label="Higher-value for restricted tilt angle search", 
                      help ='Higher-value for restricted tilt angle search \n ')             
        form.addParam('symmetry', TextParam, default='c1',
                      label='Point group symmetry',
                      help=""" See [[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry][Symmetry]]
for a description of the symmetry groups format
If no symmetry is present, give c1
""")
        form.addParam('symmetryGroupNeighbourhood', TextParam, default='', expertLevel=LEVEL_EXPERT,
                      label='Symmetry group for Neighbourhood computations',
                      help=""" If you do not know what this is leave it blank.
This symmetry will be using for compute neighboring points,
but not for sampling or reconstruction
See [[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry][Symmetry]]
for a description of the symmetry groups format
If no symmetry is present, give c1
"""
)
        form.addParam('onlyWinner', NumericListParam, default='0', 
                      label='compute only closest neighbor', expertLevel=LEVEL_EXPERT,
                      condition="symmetryGroupNeighbourhood != ''",
                      help="""This option is only relevant if symmetryGroupNeighbourhood !=''
If set to 1 only one neighbor will be computed per sampling point
You may specify this option for each iteration. 
This can be done by a sequence of numbers (for instance, "1 1 0" 
specifies 3 iterations, the first two set the value to 1 
and the last to 0. An alternative compact notation 
is ("2x1 0", i.e.,
2 iterations with value 1, and 1 with value 0).
*Note:* if there are less values than iterations the last value is reused
*Note:* if there are more values than iterations the extra value are ignored
""")     

        form.addParam('discardImages', EnumParam, 
                      choices=['None', 'maxCC', 'percentage', 'classPercentage'],
                      default=xmipp3.SELECT_NONE, display=EnumParam.DISPLAY_COMBO,
                      label='Discard images?', 
                      help=""" 
None : No images will be discarded.
maxCC  : Minimum Cross Correlation, discard images with CC below a fixed value.
percentage : Discard percentage of images with less CC.
classPercentage: Discard percentage of images in each projection direction with less CC.
Value of each option is set below.
""")
        form.addParam('minimumCrossCorrelation', NumericListParam, default='0.1', 
                      label='discard image if CC below', 
                      condition='discardImages==%d' % xmipp3.SELECT_MAXCC,
                      help=""" 
Discard images with cross-correlation (CC) below this value.
Provide a sequence of numbers (for instance, "0.3 0.3 0.5 0.5" specifies 4 iterations,
the first two set the value to 0.3, then two with 0.5.
An alternative compact notation would be ("2x0.3 2x0.5").
*Note:* if there are less values than iterations the last value is reused
*Note:* if there are more values than iterations the extra value are ignored
""")
         
        form.addParam('discardPercentage', NumericListParam, default='10', 
                      label='discard image percentage with less CC',
                      condition='discardImages==%d'%xmipp3.SELECT_PERCENTAGE,
                      help=""" 
Discard this percentage of images with less cross-correlation (CC)
Provide a sequence of numbers (for instance, "20 20 10 10" specifies 4 iterations,
the first two set the value to 20%, then two with 10%
An alternative compact notation would be ("2x20 2x10").
*Note:* if there are less values than iterations the last value is reused
*Note:* if there are more values than iterations the extra value are ignored
Set to zero to prevent discarding any images
""")
         
        form.addParam('discardPercentagePerClass', NumericListParam, default='10', 
                      label='discard image percentage in class with less CC',
                      condition='discardImages==%d'%xmipp3.SELECT_CLASSPERCENTAGE,
                      help=""" 
Discard this percentage of images in each class(projection direction)
with less cross-correlation (CC)    
Provide a sequence of numbers (for instance, "20 20 10 10" specifies 4 iterations,
the first two set the value to 20%, then two with 10%
An alternative compact notation would be ("2x20 2x10").
*Note:* if there are less values than iterations the last value is reused
*Note:* if there are more values than iterations the extra value are ignored
Set to zero to prevent discarding any images
""")
     
        form.addParam('doScale', BooleanParam, default=False,
                      label="Perform scale search?",  expertLevel=LEVEL_EXPERT,
                      help=' If true perform scale refinement. (UNDER DEVELOPMENT!!!!) \n  ')

        form.addParam('scaleStep', NumericListParam, default=1, condition='doScale',
                      label='Step scale factors size',
                      help='''Scale step factor size (1 means 0.01 in/de-crements arround 1).
Provide a sequence of numbers (for instance, "1 1 .5 .5" specifies 4 iterations,
the first two set the value to 1%, then two with .5%
An alternative compact notation would be ("2x1 2x0.5").
*Note:* if there are less values than iterations the last value is reused
*Note:* if there are more values than iterations the extra value are ignored
Set to zero to prevent discarding any images''')  

        form.addParam('ScaleNumberOfSteps', NumericListParam, default=3, condition='doScale',
                      label='Number of scale steps',
                      help=""" 
Number of scale steps.
With default values (ScaleStep='1' and ScaleNumberOfSteps='3'): 1 +/-0.01 | +/-0.02 | +/-0.03.    
With values ScaleStep='2' and ScaleNumberOfSteps='4' it performs a scale search over:
 1 +/-0.02 | +/-0.04 | +/-0.06 | +/-0.08.    
In general scale correction should only be applied to the last iteration. Do not use it unless
your data is fairly well aligned.
""")  

        form.addParam('projMatchingExtra', StringParam, default='',
                      label='Additional options for Projection_Matching', expertLevel=LEVEL_EXPERT,
                      help=""" For details see:
[[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Projection_matching][projection matching]] and
[[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Mpi_projection_matching][mpi projection matching]]
try -Ri xx -Ro yy for restricting angular search (xx and yy are
the particle inner and outter radius)
""")
        #DoSaveImagesAssignedToClasses    you can get this information in visualize
        form.addSection(label='2D re-alignment of classes', expertLevel=LEVEL_EXPERT)
        
        form.addParam('performAlign2D', BooleanParam, default=False,
                      label='Perform 2D re-alignment', expertLevel=LEVEL_EXPERT)
        
        form.addParam('doAlign2D', NumericListParam, default='0', condition='performAlign2D',
                      label='Perform 2D re-alignment of classes?', expertLevel=LEVEL_EXPERT,
                      help=""" After performing a 3D projection matching iteration, each of the
subsets of images assigned to one of the library projections is
re-aligned using a 2D-alignment protocol.
This may serve to remove model bias.
For details see:
[[http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Align2d][align 2d]]
Note that you cannot combine this option with CTF-correction!
You may specify this option for each iteration. 
This can be done by a sequence of 0 or 1 numbers (for instance, "1 1 0 0" 
specifies 4 iterations, the first two applied alig2d while the last 2
dont. an alternative compact notation is 
is ("2x1 2x0", i.e.,
2 iterations with value 1, and 2 with value 0).
*Note:*if there are less values than iterations the last value is reused
*Note:* if there are more values than iterations the extra value are ignored
*IMPORTANT:* if you set this variable to 0 the output  of the projection
muching step will be copied as output of align2d
""")
        
        form.addParam('align2DIterNr', NumericListParam, default='4', condition='performAlign2D',
                      label='Number of align2d iterations:', expertLevel=LEVEL_EXPERT,
                      help=""" Use at least 3 iterations
The number of align iteration may change in each projection matching iteration
Ffor instance, "4 4 3 3 " 
specifies 4 alig2d iterations in the first projection matching iteration 
and  two 3 alig2d iteration in the last 2 projection matching iterations.
 An alternative compact notation 
is ("2x4 2x3", i.e.,
2 iterations with value 4, and 2 with value 3).
*Note:* if there are less values than iterations the last value is reused
*Note:* if there are more values than iterations the extra value are ignored
""")        
  
  
        form.addParam('align2dMaxChangeOffset', NumericListParam, default='2x1000 2x10', 
                      condition='performAlign2D',
                      label='Maximum change in origin offset (+/- pixels)', expertLevel=LEVEL_EXPERT,
                      help="""Maximum change in shift  (+/- pixels)
You must specify this option for each iteration. 
This can be done by a sequence of numbers (for instance, "1000 1000 10 10 " 
specifies 4 iterations, the first two set the value to 1000 (no restriction)
and the last two to 10degrees. An alternative compact notation 
is ("2x1000 2x10", i.e.,
2 iterations with value 1000, and 2 with value 10).
*Note:* if there are less values than iterations the last value is reused
*Note:* if there are more values than iterations the extra value are ignored
""")    
  
        form.addParam('align2dMaxChangeRot', NumericListParam, default='2x1000 2x20', 
                      condition='performAlign2D',
                      label='Maximum change in rotation (+/- degrees)', expertLevel=LEVEL_EXPERT,
                      help="""Maximum change in shift  (+/- pixels)
You must specify this option for each iteration. 
This can be done by a sequence of numbers (for instance, "1000 1000 10 10 " 
specifies 4 iterations, the first two set the value to 1000 (no restriction)
and the last two to 10degrees. An alternative compact notation 
is ("2x1000 2x10", i.e.,
2 iterations with value 1000, and 2 with value 10).
*Note:* if there are less values than iterations the last value is reused
*Note:* if there are more values than iterations the extra value are ignored
""")     
        
        
        form.addSection(label='3D Reconstruction')
        
        form.addParam('reconstructionMethod', EnumParam, expertLevel=LEVEL_EXPERT,
                      choices=['fourier', 'art', 'wbp'],
                      default=xmipp3.RECONSTRUCT_FOURIER, display=EnumParam.DISPLAY_COMBO,
                      label='Reconstruction method', 
                      help=""" Select what reconstruction method to use.
fourier: Fourier space interpolation (with griding).
art: Agebraic reconstruction technique
wbp : Weight back project method.
""")
        
        form.addParam('fourierMaxFrequencyOfInterest', DigFreqParam, default=0.25,
                      condition='reconstructionMethod == %d' % xmipp3.RECONSTRUCT_FOURIER,
                      label='Initial maximum frequency', expertLevel=LEVEL_EXPERT,
                      help=""" This number is only used in the first iteration. 
From then on, it will be set to resolution computed in the resolution section
""")         
        form.addParam('fourierReconstructionExtraCommand', StringParam, default='',
                      condition='reconstructionMethod == %d' % xmipp3.RECONSTRUCT_FOURIER,
                      label='Additional parameters for fourier', expertLevel=LEVEL_EXPERT,
                      help=""" For details see:
[[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Fourier][fourier]]
""")          
        
        form.addParam('artLambda', NumericListParam, default='0.2', 
                      condition='reconstructionMethod == %d' % xmipp3.RECONSTRUCT_ART,
                      label='Values of lambda for ART', expertLevel=LEVEL_EXPERT,
                      help=""" *IMPORTANT:* ou must specify a value of lambda for each iteration even
if ART has not been selected.
*IMPORTANT:* NOte that we are using the WLS version of ART that 
uses geater lambdas than the plain art.
See for details:
[[http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Art][xmipp art]]
You must specify this option for each iteration. 
This can be done by a sequence of numbers (for instance, ".1 .1 .3 .3" 
specifies 4 iterations, the first two set the value to 0.1 
(no restriction)
and the last  two to .3. An alternative compact notation 
is ("2x.1 2x.3").
*Note:* if there are less values than iterations the last value is reused
*Note:* if there are more values than iterations the extra value are ignored
""")   
        
        form.addParam('artReconstructionExtraCommand', StringParam, default='-k 0.5 -n 10 ',
                      condition='reconstructionMethod == %d' % xmipp3.RECONSTRUCT_ART,
                      label='Additional parameters for ART', expertLevel=LEVEL_EXPERT,
                      help=""" For details see:
[[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Art][xmipp art]]
""")          
        
        form.addParam('wbpReconstructionExtraCommand', StringParam, default='',
                      condition='reconstructionMethod == %d' % xmipp3.RECONSTRUCT_WBP,
                      label='Additional parameters for WBP', expertLevel=LEVEL_EXPERT,
                      help=""" For details see:
[[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Wbp][xmipp wbp]]
""")                  

        form.addParam('doComputeResolution', NumericListParam, default='1',
                      label='Compute resolution?', expertLevel=LEVEL_EXPERT,
                      help=""" For details see:
    [[http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Resolution][xmipp resolution]]
    Set to 1 to compute resolution and to 0 if you do not want to compute it.
""")
    
        form.addParam('doSplitReferenceImages', NumericListParam, default='1',
                      label='Split references averages?', expertLevel=LEVEL_EXPERT,
                      condition="doComputeResolution!='0'",
                      help="""In theory each reference average should be splited
in two when computing the resolution. In this way each
projection direction will be represented in each of the
subvolumes used to compute the resolution. A much faster
but less accurate approach is to split the 
proyection directions in two but not the averages. We
recomend the first approach for small volumes and the second for
large volumes (especially when using small angular
sampling rates.
*IMPORTANT:* the second option has ONLY been implemented for FOURIER
reconstruction method. Other reconstruction methods require this
flag to be set to True
 You may specify this option for each iteration. 
 This can be done by a sequence of 0 or 1 numbers (for instance, "1 1 0 0" 
 specifies 4 iterations, the first two split the images   while the last 2
 don't. an alternative compact notation is 
 is ("2x1 2x0", i.e.,
 2 iterations with value 1, and 2 with value 0).
 *Note:* if there are less values than iterations the last value is reused
 *Note:* if there are more vapplications/scripts/protocols/new_protocol_projmatch.pyalues than iterations the extra value are ignored
""")            

        form.addParam('doLowPassFilter', BooleanParam, default=True,
                      label="Low-pass filter the reference?",  expertLevel=LEVEL_ADVANCED,
                      help="""If set to true, the volume will be filtered at a frecuency equal to
the  resolution computed with a FSC=0.5 threshold, possibly 
plus a constant provided by the user in the next input box. 

If set to false, then the filtration will be made at the constant 
value provided by the user in the next box (in digital frequency, 
i.e. pixel-1: minimum 0, maximum 0.5) 
""")

        form.addParam('constantToAddToFiltration', NumericListParam, default='0.1 0',
                      label='Constant to be added to the estimated resolution', expertLevel=LEVEL_ADVANCED,
                      condition="doLowPassFilter!='0'",
                      help=""" The meaning of this field depends on the previous flag.
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
*Note:* if there are less values than iterations the last value is reused
*Note:* if there are more values than iterations the extra value are ignored
""")
        form.addParam('ConstantToAddToMaxReconstructionFrequency', NumericListParam, default='0.1',
                      label='Constant to be added to the reconstruction maximum frequency', expertLevel=LEVEL_ADVANCED,
                      condition="doLowPassFilter!='0'",
                      help=""" The meaning of this field depends on the UseFscForFilter flag.
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
*Note:* if there are less values than iterations the last value is reused
*Note:* if there are more values than iterations the extra value are ignored
""")
        
#        
#        form.addSection(label='ML3D copy')       
#        
#        form.addParam('numberOfSeedsPerRef', IntParam, default=1,
#                      label='Number of seeds per reference',
#                      help='The total number of seeds generated will be the number of provided '
#                      'reference volumes times the number of seeds per reference. '
#                      'If you provide 2 initial volumes and 3 seeds per referece you will '
#                      'produce 6 3D maps.')
#        form.addParam('angularSampling', IntParam, default=10,
#                      label='Angular sampling for classification',
#                      help='Fine samplings take huge amounts of CPU and memory. '
#                      'Therefore, in general, dont use samplings finer than 10 degrees.')    
#        form.addParam('symmetry', TextParam, default='c1',
#                      label='Point group symmetry',
#                      help='Number of ML(F)3D iterations to perform.')        
#        form.addParam('restartIter', IntParam, default=0,
#                      expertLevel=LEVEL_ADVANCED,
#                      label='Restart after iteration',
#                      help='For previous runs that stopped before convergence, '
#                      'resume the calculations after the completely finished iteration, '
#                      'i.e. including all 3D reconstructions. '
#                      'Note that all flags about grey-scale correction, filtering and '
#                      'seed generation will be ignored if a value larger than 0 is given, '
#                      'since this option only concerns the ProjMatch classification part.')   
#        form.addParam('extraParams', TextParam,
#                      expertLevel=LEVEL_ADVANCED,
#                      label='Additional parameters',
#                      help='Additional xmipp_ml(f)_refine3d parameters.')                  
#        form.addSection(label='MLF parameters', questionParam='doMlf')        
#        form.addParam('doMlf', BooleanParam, default=False,
#                      label='Use MLF2D instead of ML2D')
#        form.addParam('doCorrectAmplitudes', BooleanParam, default=True,
#                      label='Use CTF-amplitude correction inside MLF?',
#                      help='If set to <Yes>, the input images file should contain '
#                           'the CTF information for each image.')
#        form.addParam('highResLimit', IntParam, default=20,
#                      label='High-resolution limit (in Angstroms)',
#                      help='No frequencies higher than this limit will be taken into account. '
#                      'If zero is given, no limit is imposed.')     
#        form.addParam('areImagesPhaseFlipped', BooleanParam, default=True,
#                      label='Are the images CTF phase flipped?',
#                      help='You can run MLF with or without having phase flipped the images.')     
#        form.addParam('initialMapIsAmplitudeCorrected', BooleanParam, default=False,
#                      label='Are initial references CTF-amplitude corrected?',
#                      help='If coming from programs other than xmipp_mlf_refine3d this is '
#                      'usually not the case. If you will perform a grey-scale correction, '
#                      'this parameter becomes irrelevant as the output maps never have the '
#                      'CTF-amplitudes corrected.')     
#        form.addParam('seedsAreAmplitudeCorrected', BooleanParam, default=False,
#                      expertLevel=LEVEL_ADVANCED,
#                      label='Are the seeds CTF-amplitude corrected?',
#                      help='This option is only relevant if you provide your own seeds! '
#                      'If the seeds are generated automatically, this parameter becomes '
#                      'irrelevant as they will always be amplitude-corrected.') 
#        form.addSection(label='3D Reconstruction', expertLevel=LEVEL_ADVANCED)    
#        form.addParam('reconstructionMethod', EnumParam, choices=['fourier', 'wlsART'], 
#                      default=0, label='Reconstruction method', display=EnumParam.DISPLAY_LIST,
#                      expertLevel=LEVEL_ADVANCED, 
#                      help='Choose between wslART or fourier.')
#        form.addParam('aRTExtraParams', TextParam,
#                      condition='reconstructionMethod==1', expertLevel=LEVEL_ADVANCED,
#                      label='Extra parameters',
#                      help='Additional reconstruction parameters for ART.')  
#        form.addParam('fourierExtraParams', TextParam, 
#                      condition='reconstructionMethod==0', expertLevel=LEVEL_ADVANCED,
#                      label='Extra parameters',
#                      help='The Fourier-interpolation reconstruction method is much faster than wlsART '
#                      'and may give similar results. It however is not guaranteed to optimize the '
#                      'likelihood function. This is an experimental feature. One may limit the '
#                      'maximum resolution of the fourier-interpolation using -max_resolution 0.3 '
#                      '(to 0.3 digital frequency). Use the extra parameter entry below for that.')  
        
        form.addParallelSection(threads=1, mpi=8)
        
    def getProgramId(self):
        progId = "ml"
        if self.doMlf:
            progId += "f" 
        return progId   
         
    def _defineSteps(self):

        self.ParamsDict = {}
        self.ParamsDict['ProgId'] = self.getProgramId()
        
        #TODO: Retrieve initial volume path from object (now is a text) and convert if needed
        refMd = self.ini3DrefVolumes.get()
        initVols = self.ParamsDict['InitialVols'] = self._getExtraPath('initial_volumes.stk')
        self.mdVols = xmipp.MetaData(refMd)
        
        # Convert input images if necessary
        self.inputImgs = self.inputImages.get()        
        imgsFn = self._insertConvertStep('inputImgs', XmippSetOfImages,
                                         self._getPath('input_images.xmd'))
        self.imgMd = self.ParamsDict['ImgMd'] = imgsFn
        
        #Get sampling rate from input images
        self.samplingRate = self.inputImages.get().getSamplingRate()
        
        self._insertFunctionStep('copyVolumes', refMd, initVols)
        
        if self.doCorrectGreyScale.get():
            self.insertCorrectGreyScaleSteps()
                    
        if self.doLowPassFilter.get():
            self.insertFilterStep()
            
        if self.numberOfSeedsPerRef.get() > 1:
            self.insertGenerateRefSteps()
            
        self.insertProjMatchStep(self.imgMd, self.workingDir.get() + '/', self.ParamsDict['InitialVols'], 
                            self.numberOfIterations.get(), self.seedsAreAmplitudeCorrected.get())
        
        self._insertFunctionStep('renameOutput', self.workingDir.get(), self.getProgramId())
                
        self._insertFunctionStep('createOutput')

    def copyVolumes(self, inputMd, outputStack):
        ''' This function will copy input references into a stack in working directory'''
        if exists(outputStack):
            os.remove(outputStack)
        md = xmipp.MetaData(inputMd)
        img = xmipp.Image()
        for i, idx in enumerate(md):
            img.read(md.getValue(xmipp.MDL_IMAGE, idx))
            img.write('%d@%s' % (i + 1, outputStack))

        
    def insertCorrectGreyScaleSteps(self):
        ''' Correct the initial reference greyscale '''
        cgsDir = self._getPath('CorrectGreyscale')
        makePath(cgsDir)
        volStack = self.ParamsDict['InitialVols'] = self._getExtraPath('corrected_volumes.stk')
        # Grey-scale correction always leads to an amplitude uncorrected map
        self.initialMapIsAmplitudeCorrected.set(False)
        index = 1
        outputVol = ''
        for idx in self.mdVols:
            volDir = join(cgsDir, 'vol%03d' % index)
            projs = join(volDir, 'projections')
            makePath(volDir)
            outputVol = "%(index)d@%(volStack)s" % locals()
            corrRefsRoot = join(volDir, 'corrected_refs')
            self.ParamsDict.update({
                'inputVol': self.mdVols.getValue(xmipp.MDL_IMAGE, idx),
                'outputVol': outputVol,
                'projRefs': projs + ".stk",
                'docRefs': projs + ".doc",
                'corrRefsRoot':corrRefsRoot ,
                'corrRefs': corrRefsRoot + '_Ref3D_001.stk',
                'projMatch': join(volDir, "proj_match.doc"),
                'projMatchSampling': self.projMatchSampling.get(),
                'symmetry': self.symmetry.get(),
                'numberOfThreads': self.numberOfThreads.get()
                })
            self.mdVols.setValue(xmipp.MDL_IMAGE, outputVol, idx)
            self.ParamsStr = ' -i %(inputVol)s --experimental_images %(ImgMd)s -o %(projRefs)s' + \
                    ' --sampling_rate %(projMatchSampling)f --sym %(symmetry)s' + \
                    'h --compute_neighbors --angular_distance -1' 
                       
            self._insertRunJobStep('xmipp_angular_project_library', self.ParamsStr % self.ParamsDict)

            self.ParamsStr = '-i %(ImgMd)s -o %(projMatch)s --ref %(projRefs)s' 
            self._insertRunJobStep('xmipp_angular_projection_matching', self.ParamsStr % self.ParamsDict)
 
            self.ParamsStr = '-i %(projMatch)s --lib %(docRefs)s -o %(corrRefsRoot)s'
            self._insertRunJobStep('xmipp_angular_class_average', self.ParamsStr % self.ParamsDict)

            self.ParamsStr = '-i %(projMatch)s -o %(outputVol)s --sym %(symmetry)s --weight --thr %(numberOfThreads)d'
            self._insertRunJobStep('xmipp_reconstruct_fourier', self.ParamsStr % self.ParamsDict)
            index += 1

    def insertFilterStep(self):
        volStack = self.ParamsDict['FilteredVols'] = self._getExtraPath('filtered_volumes.stk')
        index = 1
        outputVol = ''
        for idx in self.mdVols:
            outputVol = "%(index)d@%(volStack)s" % locals()
            self.mdVols.setValue(xmipp.MDL_IMAGE, outputVol, idx)
            index += 1
        self.ParamsStr = '-i %(InitialVols)s -o %(FilteredVols)s --fourier low_pass %(lowPassFilter)f --sampling %(samplingRate)f'
        self.ParamsDict.update({
                                'lowPassFilter':self.lowPassFilter.get(),
                                'samplingRate':self.samplingRate
                                })
        self._insertRunJobStep('xmipp_transform_filter', self.ParamsStr % self.ParamsDict)
        self.ParamsDict['InitialVols'] = self.ParamsDict['FilteredVols']

    def insertGenerateRefSteps(self):
        ''' Generate more reference volumes than provided in input reference '''
        grDir = self._getPath('GeneratedReferences')
        # Create dir for seeds generation
        makePath(grDir)
        # Split images metadata
        nvols = self.ParamsDict['NumberOfVols'] = self.mdVols.size() * self.numberOfSeedsPerRef.get()
        sroot = self.ParamsDict['SplitRoot'] = join(grDir, 'images')
        self.ParamsStr = '-i %(ImgMd)s -n %(NumberOfVols)d --oroot %(SplitRoot)s'
        files = ['%s%06d.xmd' % (sroot, i) for i in range(1, nvols+1)]        
        self._insertRunJobStep('xmipp_metadata_split', self.ParamsStr % self.ParamsDict, numberOfMpi=1)
        
        volStack = self.ParamsDict['InitialVols'] = self._getExtraPath('generated_volumes.stk') 
        index = 1
        copyVols = []
        for idx in self.mdVols:
            for i in range(self.numberOfSeedsPerRef.get()):
                outputVol = "%d@%s" % (index, volStack)
                generatedVol = join(grDir, "vol%03dextra/iter%03d/vol%06d.vol" % (index, 1, 1))
                copyVols.append((outputVol, generatedVol))
                self.insertProjMatchStep(files[index-1], join(grDir, 'vol%03d' % index), self.mdVols.getValue(xmipp.MDL_IMAGE, idx), 1, 
                                    self.initialMapIsAmplitudeCorrected)
                #self.mdVols.setValue(MDL_IMAGE, outputVol, idx)
                index += 1
                
        for outVol, genVol in copyVols:
            self.ParamsDict.update({'outVol': outVol, 'genVol':genVol})
            self.ParamsStr = '-i %(genVol)s -o %(outVol)s'
            self._insertRunJobStep('xmipp_image_convert', self.ParamsStr % self.ParamsDict, numberOfMpi=1)
            
        # Seed generation with MLF always does amplitude correction
        self.seedsAreAmplitudeCorrected.set(True)

    def insertProjMatchStep(self, inputImg, oRoot, initialVols, numberOfIters, amplitudCorrected):
        self.ParamsDict.update({
                         '_ImgMd': inputImg,
                         '_ORoot': oRoot,
                         '_InitialVols': initialVols,
                         '_NumberOfIterations': numberOfIters,
                         'symmetry':self.symmetry.get(),
                         'angularSampling': self.angularSampling.get(),
                         'extraParams': self.extraParams.get(),
                         'numberOfThreads': self.numberOfThreads.get(),
                         'highResLimit': self.highResLimit.get(),
                         'samplingRate': self.samplingRate,
                         'reconstructionMethod': self.getEnumText('reconstructionMethod'),
                         'aRTExtraParams': self.aRTExtraParams.get(),
                         'fourierExtraParams': self.fourierExtraParams.get()
                        })
        self.ParamsStr = "-i %(_ImgMd)s --oroot %(_ORoot)s --ref %(_InitialVols)s --iter %(_NumberOfIterations)d " + \
                         "--sym %(symmetry)s --ang %(angularSampling)s"
        if self.aRTExtraParams.hasValue(): self.ParamsStr += " %(extraParams)s"
#        if self.NumberOfReferences > 1:
#            self.ParamsStr += " --nref %(NumberOfReferences)s"
        if self.numberOfThreads.get() > 1:
            self.ParamsStr += " --thr %(numberOfThreads)d"
        if self.doNorm.get():
            self.ParamsStr += " --norm"
        
        if self.doMlf.get():
            if not self.doCorrectAmplitudes.get():
                self.ParamsStr += " --no_ctf"
            if not self.areImagesPhaseFlipped.get():
                self.ParamsStr += " --not_phase_flipped"
            if not amplitudCorrected:
                self.ParamsStr += " --ctf_affected_refs"
            if self.highResLimit > 0:
                self.ParamsStr += " --limit_resolution 0 %(highResLimit)f"
            self.ParamsStr += ' --sampling_rate %(samplingRate)f'

        self.ParamsStr += " --recons %(reconstructionMethod)s "
        
        if self.reconstructionMethod.get() == self.WLSART:
            if self.aRTExtraParams.hasValue(): self.ParamsStr += " %(aRTExtraParams)s"
        else:
            if self.fourierExtraParams.hasValue(): self.ParamsStr += " %(fourierExtraParams)s" 
            
        self._insertRunJobStep('xmipp_%s_refine3d' % self.getProgramId(), self.ParamsStr % self.ParamsDict)
        
    def renameOutput(self, WorkingDir, ProgId):
        ''' Remove ml2d prefix from:
            ml2dclasses.stk, ml2dclasses.xmd and ml2dimages.xmd'''
        prefix = '%s2d' % ProgId
        for f in ['%sclasses.stk', '%sclasses.xmd', '%simages.xmd']:
            f = join(WorkingDir, f % prefix)
            nf = f.replace(prefix, '')
            shutil.move(f, nf)
                                                    
    def createOutput(self):
        lastIter = 'iter%03d' % self.numberOfIterations.get()
        volumes = XmippSetOfVolumes(self._getExtraPath(lastIter + '/' + 'iter_volumes.xmd'))
        self._defineOutputs(outputVolumes=volumes)

    def _summary(self):
        summary = []
        summary.append("Input images:  %s" % self.inputParticles.get().getNameId())
        if self.doCTFCorrection.get():
            suffix = "with CTF correction "
        else:
            suffix = "ignoring CTF effects "
        summary.append("Using a ML in *Fourier-space* " + suffix)
         
        summary.append("Reference volumes(s): [%s]" % self.input3DReferences.get())

#        if self.numberOfSeedsPerRef.get() > 1:
#            summary.append("Number of references per volume: <%d>" % self.numberOfSeedsPerRef.get())
           
        # TODO: Add information at info@iter_classes.xmd from last iteration
        
        if not hasattr(self, 'outputVolumes'):
            summary.append("Output volumes not ready yet.")
        else:
            summary.append("Output volumes: %s" % self.outputVolumes.get())
        return summary
    
    def _validate(self):
        validateMsgs = []
#        # Volume references cannot be empty (REMOVE this when we use a pointer to a set instead of a path)
#        if not self.ini3DrefVolumes.hasValue():
#            validateMsgs.append('Please provide an initial reference volume(s).')
#        # Correct grey scale needs mpi to b e run
#        if self.doCorrectGreyScale.get() and self.numberOfMpi < 2:
#            validateMsgs.append('Correct grey scale needs mpi to be run.')
            
        #TODO: Check images dimension when it is implemented on SetOFImages class
        return validateMsgs