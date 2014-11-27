# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *              Josue Gomez Blanco (jgomez@cnb.csic.es)
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
This sub-package contains wrapper around XmippProtPreprocessVolumes protocol
"""
from pyworkflow.em import *  
from pyworkflow.utils import *  
import xmipp
from protocol_process import XmippProcessParticles, XmippProcessVolumes
from pyworkflow.utils.path import cleanPath
from pyworkflow.em.constants import *
from constants import *
from convert import locationToXmipp



class XmippPreprocessHelper():
    """ 
    Helper class that contains some Protocol utilities methods
    used by both  XmippProtPreprocessParticles and XmippProtPreprocessVolumes.
    """
    
    #--------------------------- DEFINE param functions --------------------------------------------
    @classmethod
    def _defineProcessParams(cls, form):
        # Invert Contrast
        form.addParam('doInvert', BooleanParam, default=False,
                      label='Invert contrast',
                      help='Invert the contrast if your particles are black over a white background.')
        # Threshold
        form.addParam('doThreshold', BooleanParam, default=False,
                      label="Threshold",
                      help='Remove voxels below a certain value.')
        form.addParam('thresholdType', EnumParam, condition='doThreshold',
                      choices=['abs_below', 'below', 'above'],
                      default=MASK_FILL_VALUE,
                      label="Fill with ", display=EnumParam.DISPLAY_COMBO,
                      help='Select how are you going to fill the pixel values outside the mask. ')
        form.addParam('threshold', FloatParam, default=0,
                      label="Threshold value", condition='doThreshold',
                      help='Grey value below which all voxels should be set to 0.')
        form.addParam('fillType', EnumParam, condition='doThreshold',
                      choices=['value', 'binarize', 'avg'],
                      default=FILL_VALUE,
                      label="Substitute by", display=EnumParam.DISPLAY_COMBO,
                      help='If you select: value: Selected are substitute by a desired value.\n'
                           '            binarize: Selected are set to 0, non-selected to 1.\n'
                           '                 avg: Average of non-selected.')
        form.addParam('fillValue', IntParam, default=0, condition='doThreshold and fillType == %d'  % FILL_VALUE,
                      label='Fill value',
                      help=' Substitute selected pixels by this value.')

    #--------------------------- INSERT steps functions --------------------------------------------
    @classmethod
    def _insertCommonSteps(cls, protocol, changeInserts):
        if protocol.doInvert:
            args = protocol._argsInvert()
            if protocol.isFirstStep:
                protocol.isFirstStep = False
            protocol._insertFunctionStep("invertStep", args, changeInserts)

        if protocol.doThreshold:
            args = protocol._argsThreshold()
            if protocol.isFirstStep:
                protocol.isFirstStep = False
            protocol._insertFunctionStep("thresholdStep", args, changeInserts)

    #--------------------------- UTILS functions ---------------------------------------------------
    @classmethod
    def _argsCommonInvert(cls):
        args = ' --mult -1'
        return args
    
    @classmethod
    def _argsCommonThreshold(cls, protocol):
        args = " --select below %f" % protocol.threshold.get()
        fillStr = protocol.getEnumText('fillType')
        args += " --substitute %s " % fillStr
        
        if protocol.fillType == MASK_FILL_VALUE:
            args += " %f" % protocol.fillValue.get()
        return args


class XmippProtPreprocessParticles(XmippProcessParticles,XmippPreprocessHelper):
    """ Preprocess a set of particles. You can remove dust, normalize, apply threshold, etc """
    _label = 'preprocess particles'

    # Automatic Particle rejection enum
    REJ_NONE = 0
    REJ_MAXZSCORE = 1
    REJ_PERCENTAGE =2
    
    def __init__(self, **args):
        XmippProcessParticles.__init__(self)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        form.addParam('doRemoveDust', BooleanParam, default=False,
                      label='Dust removal', 
                      help='Sets pixels with unusually large values to random values from a Gaussian '
                      'with zero-mean and unity-standard deviation.')
        form.addParam('thresholdDust', FloatParam, default=3.5, condition='doRemoveDust',
                      label='Threshold for dust removal',
                      help='Pixels with a signal higher or lower than this value times the standard '
                      'deviation of the image will be affected. For cryo, 3.5 is a good value.'
                      'For high-contrast negative stain, the signal itself may be affected so '
                      'that a higher value may be preferable.',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('doNormalize', BooleanParam, default=False,
                      label='Normalize', 
                      help='It subtract a ramp in the gray values and normalizes so that in the '
                      'background there is 0 mean and standard deviation 1.')
        form.addParam('normType', EnumParam, choices=['OldXmipp','NewXmipp','Ramp'], 
                      default=2, condition='doNormalize', display=EnumParam.DISPLAY_COMBO,
                      label='Normalization type', 
                      help='OldXmipp (mean(Image)=0, stddev(Image)=1). \n  '
                           'NewXmipp (mean(background)=0, stddev(background)=1) \n  '
                           'Ramp (subtract background+NewXmipp). \n',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('backRadius', IntParam, default=-1, condition='doNormalize',
                      label='Background radius',
                      help='Pixels outside this circle are assumed to be noise and their stddev '
                      'is set to 1. Radius for background circle definition (in pix.). '
                      'If this value is 0, then half the box size is used.', 
                      expertLevel=LEVEL_ADVANCED)
        XmippPreprocessHelper._defineProcessParams(form)
#         form.addParam('autoParRejection', EnumParam, choices=['None', 'MaxZscore', 'Percentage'],
#                       label="Automatic particle rejection", default=REJ_NONE,
#                       display=EnumParam.DISPLAY_COMBO, expertLevel=LEVEL_EXPERT,
#                       help='How to automatically reject particles. It can be none (no rejection), '
#                       'maxZscore (reject a particle if its Zscore is larger than this value), '
#                       'Percentage (reject a given percentage in each one of the screening criteria). ')
#         form.addParam('maxZscore', IntParam, default=3, condition='autoParRejection==1',
#                       label='Maximum Zscore', expertLevel=LEVEL_EXPERT,
#                       help='Maximum Zscore.', validators=[Positive])      
#         form.addParam('percentage', IntParam, default=5, condition='autoParRejection==2',
#                       label='Percentage (%)', expertLevel=LEVEL_EXPERT,
#                       help='Percentage.', validators=[Range(0, 100, error="Percentage must be between 0 and 100.")])        
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertProcessStep(self):
        self.isFirstStep = True
        # this is for when the options selected has changed and the protocol is resumed
        changeInserts = [self.doRemoveDust, self.doNormalize, self.doInvert, self.doThreshold]
        
        if self.doRemoveDust:
            args = self._argsRemoveDust()
            if self.isFirstStep:
                self.isFirstStep = False
            self._insertFunctionStep("removeDustStep", args, changeInserts)
        
        if self.doNormalize:
            args = self._argsNormalize()
            if self.isFirstStep:
                self.isFirstStep = False
            self._insertFunctionStep("normalizeStep", args, changeInserts)
        
        XmippPreprocessHelper._insertCommonSteps(self, changeInserts)
        
#         if self.getEnumText('autoParRejection') != 'None':
#             self._insertFunctionStep("rejectionStep", outputFn, outputMd)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def removeDustStep(self, args, changeInserts):
        self.runJob('xmipp_transform_filter', args)
    
    def normalizeStep(self, args, changeInserts):
        self.runJob("xmipp_transform_normalize", args % locals())
    
    def sortImages(self, outputFn, outputMd):
        pass
#         from xmipp import MetaDataInfo
#         #(Xdim, Ydim, Zdim, Ndim, _) = MetaDataInfo(inputFile)
#         args=""
#         # copy file to run path
#         self.outputMd = String(self._getPath(replaceBaseExt(inputFile, 'xmd')))
#         self.outputMd._objDoStore = True
#         if inputFile != self.outputMd.get():
#             copyFile(inputFile, self.outputMd.get())
#         if self.autoParRejection.get()==REJ_MAXZSCORE:
#             args+=" --zcut "+str(self.maxZscore.get())
#         elif self.autoParRejection.get()==REJ_PERCENTAGE:
#             args+=" --percent "+str(self.percentage.get())
#         #if Ndim > 0:
#         self.runJob("xmipp_image_sort_by_statistics", "-i " + self.outputMd.get() + " --addToInput"+args)
#     
#     def convertInputStep(self):
#         """ convert if necessary"""
#         imgSet = self.inputParticles.get()
#         imgSet.writeStack(self.outputStk)
#     
#     def createOutputStep(self):
#         inImgSet = self.inputParticles.get()
#         outImgSet = self._createSetOfParticles()
#         outImgSet.copyInfo(inImgSet)
#         
#         for i, img in enumerate(inImgSet):
#             j = i + 1
#             img.setLocation(j, self.outputStk)
#             outImgSet.append(img)
#         
#         self._defineOutputs(outputParticles=outImgSet)
#         self._defineTransformRelation(inImgSet, self.outputParticles)
#     
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        validateMsgs = []
        
        if self.doNormalize.get() and self.normType.get() != 0:
            size = self._getSize()
            if self.backRadius.get() > size:
                validateMsgs.append('Set a valid Background radius less than %d' % size)
        return validateMsgs
    
    def _summary(self):
        summary = []
        summary.append("Input particles: %s" % self.inputParticles.get().getFileName())
        
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        else:
            summary.append("Dust removal: %s" % self.doRemoveDust)
            summary.append("Normalize the background: %s" % self.doNormalize)
            summary.append("Invert contrast: %s" % self.doInvert)
            summary.append("Remove voxels with threshold: %s" % self.doThreshold)
        return summary
    
    def _methods(self):
        methods = []
        methods.append("Input particles: %s" % self.inputParticles.get().getFileName())
        
        if hasattr(self, 'outputParticles'):
            methods.append("A set of %d particles was preprocessed." % self.outputParticles.getSize())
            if self.doNormalize:
                methods.append("The background was normalized with %s method." % self.getEnumText('normType'))
            if self.doInvert:
                methods.append("The contrast was inverted")
            if self.doThreshold:
                methods.append("Pixels with values below %f was removed" % self.threshold.get())
        return methods
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _argsRemoveDust(self):
        if self.isFirstStep:
            args = "-i %s -o %s --save_metadata_stack %s --keep_input_columns" % (self.inputFn, self.outputStk, self.outputMd)
        else:
            args = "-i %s" % self.outputStk
        args += " --bad_pixels outliers %f" % self.thresholdDust.get()
        return args
    
    def _argsNormalize(self):
        if self.isFirstStep:
            args = "-i %s -o %s --save_metadata_stack %s --keep_input_columns" % (self.inputFn, self.outputStk, self.outputMd)
        else:
            args = "-i %s" % self.outputStk
        
        normType = self.normType.get()
        bgRadius = self.backRadius.get()
        radii = self._getSize()
        if bgRadius <= 0:
            bgRadius = int(radii)
        
        if normType == "OldXmipp":
            args += " --method OldXmipp"
        elif normType == "NewXmipp":
            args += " --method NewXmipp --background circle %d" % bgRadius
        else:
            args += " --method Ramp --background circle %d" % bgRadius
        return args
    
    def _argsInvert(self):
        if self.isFirstStep:
            args = "-i %s -o %s --save_metadata_stack %s --keep_input_columns" % (self.inputFn, self.outputStk, self.outputMd)
        else:
            args = "-i %s" % self.outputStk
        #args += XmippPreprocessHelper._argsCommonInvert()
        args += self._argsCommonInvert()
        return args
    
    def _argsThreshold(self):
        if self.isFirstStep:
            args = "-i %s -o %s --save_metadata_stack %s --keep_input_columns" % (self.inputFn, self.outputStk, self.outputMd)
        else:
            args = "-i %s" % self.outputStk
        args += self._argsCommonThreshold()
        return args
    
    def _getSize(self):
        """ get the size of SetOfParticles object"""
        Xdim = self.inputParticles.get().getDimensions()[0]
        size = int(Xdim/2)
        return size


class XmippProtPreprocessVolumes(XmippProcessVolumes):
    """ Protocol for Xmipp-based preprocess for volumes """
    _label = 'preprocess volumes'
    
    # Aggregation constants
    AGG_AVERAGE=0
    AGG_SUM=1
    
    # Segmentation type
    SEG_VOXEL=0
    SEG_AMIN=1
    SEG_DALTON=2
    SEG_AUTO=3


    def __init__(self, **args):
        XmippProcessVolumes.__init__(self)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        # Change hand
        form.addParam('doChangeHand', BooleanParam, default=False,
                      label="Change hand", 
                      help='Change hand by applying a mirror along X.')
        # Randomize the phases
        form.addParam('doRandomize', BooleanParam, default=False,
                      label="Randomize phases", 
                      help='Randomize phases beyond a certain frequency.')
        # ToDo: add wizard
        form.addParam('maxResolutionRandomize', FloatParam, default=40,
                      label="Maximum Resolution", condition='doRandomize',
                      help='Angstroms.')
        # Symmetrization
        form.addParam('doSymmetrize', BooleanParam, default=False,
                      label="Symmetrize", 
                      help='Symmetrize the input model.')
        form.addParam('symmetryGroup', TextParam, default='i1',
                      label="Symmetry group", condition='doSymmetrize',
                      help='See [[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry][Symmetry]] '
                      'for a description of the symmetry groups format.\nIf no symmetry is present, set the Symmetrize field to not.')
        form.addParam('aggregation', EnumParam, choices=['Average', 'Sum'], 
                      display=EnumParam.DISPLAY_COMBO,
                      default=0, label='Aggregation mode', condition = 'doSymmetrize',
                      help='Symmetrized volumes can be averaged or summed.')
        # Adjust gray values
        form.addParam('doAdjust', BooleanParam, default=False,
                      label="Adjust gray values", 
                      help='Adjust input gray values so that it is compatible with a set of projections.') 
        form.addParam('setOfProjections', PointerParam, pointerClass='SetOfParticles',
                      label="Set of projections", condition='doAdjust',
                      help='Set of images to which the model should conform. The set of images should have the'
                      'final pixel size and the final size of the model.')
        # Segment
        form.addParam('doSegment', BooleanParam,
                      default=False, label="Segment", 
                      help='Separate the molecule from its background.')
        form.addParam('segmentationType', EnumParam, choices=['Voxel mass', 'Aminoacid mass','Dalton mass','Automatic'],
                      default=3, display=EnumParam.DISPLAY_COMBO,
                      label="Segmentation Type", condition='doSegment',
                      help='Type of segmentation.')
        form.addParam('segmentationMass', FloatParam, default=-1,
                      label="Molecule Mass", condition='doSegment and segmentationType != 3',
                      help='In automatic segmentation, set it to -1.')
        # Normalize background
        form.addParam('doNormalize', BooleanParam, default=False,
                      label="Normalize background", 
                      help='Set background to have zero mean and standard deviation 1.')
        form.addParam('backRadius', FloatParam, default=-1,
                      label="Mask Radius", condition='doNormalize',
                      help='In pixels. Set to -1 for half of the size of the volume.')
        #XmippPreprocessHelper._defineProcessParams(form)
        self._defineProcessParams(form)

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertProcessStep(self):
        self.isFirstStep = True
        # this is for when the options selected has changed and the protocol is resumed
        changeInserts = [self.doChangeHand, self.doRandomize, self.doSymmetrize, 
                         self.doAdjust, self.doSegment, self.doInvert, self.doNormalize, 
                         self.doThreshold]

        if self.doChangeHand:
            args = self._argsChangeHand()
            if self.isFirstStep:
                self.isFirstStep = False
            self._insertFunctionStep("changeHandStep", args, changeInserts)
        
        if self.doRandomize:
            args = self._argsRandomize()
            if self.isFirstStep:
                self.isFirstStep = False
            self._insertFunctionStep("randomizeStep", args, changeInserts)
        
        if self.doSymmetrize:
            args = self._argsSymmetrize()
            if self.isFirstStep:
                self.isFirstStep = False
            self._insertFunctionStep("symmetrizeStep", args, changeInserts)
        
        if self.doAdjust:
            args = self._argsAdjust()
            self._insertFunctionStep("adjustStep", args, changeInserts)
            if self.isFirstStep:
                self.isFirstStep = False

        if self.doSegment:
            args = self._argsSegment()
            self._insertFunctionStep("segmentStep", args, changeInserts)
            if self.isFirstStep:
                self.isFirstStep = False
        
        if self.doNormalize:
            args = self._argsNormalize()
            if self.isFirstStep:
                self.isFirstStep = False
            self._insertFunctionStep("normalizeStep", args, changeInserts)
        
        #XmippPreprocessHelper._insertCommonSteps(self, changeInserts)
        self._insertCommonSteps(self, changeInserts)

    #--------------------------- STEPS functions ---------------------------------------------------
    def removeDustStep(self, args, changeInserts):
        self.runJob('xmipp_transform_filter', args)
    
    def normalizeStep(self, args, changeInserts):
        self.runJob("xmipp_transform_normalize", args % locals())
        
    def changeHandStep(self, args, changeInserts):
        self.runJob("xmipp_transform_mirror", args)
    
    def randomizeStep(self, args, changeInserts):
        self.runJob("xmipp_transform_randomize_phases", args)
    
    def symmetrizeStep(self, args, changeInserts):
        self.runJob("xmipp_transform_symmetrize", args)
    
    def adjustStep(self, args, changeInserts):
        if self._isSingleInput():
            localArgs = self._adjustLocalArgs(self.inputFn, self.outputStk, args)
            self.runJob("xmipp_transform_adjust_volume_grey_levels", localArgs)
        else:
            numberOfVols = self.inputVolumes.get().getSize()
            
            for i in range(1, numberOfVols + 1):
                inputVol = locationToXmipp(i, self.inputFn)
                outputVol = locationToXmipp(i, self.outputStk)
                localArgs = self._adjustLocalArgs(self.inputFn, self.outputStk, args)
                self.runJob("xmipp_transform_adjust_volume_grey_levels", localArgs)
    
    def segmentStep(self, args, changeInserts):
        fnMask = self._getTmpPath("mask.vol")
        if self._isSingleInput():
            localArgs = self._segmentLocalArgs(self.inputFn, fnMask, args)
            maskArgs = self._segMentMaskArgs(self.inputFn, self.outputStk, fnMask)
            self._segmentVolume(localArgs, maskArgs, fnMask)
        else:
            numberOfVols = self.inputVolumes.get().getSize()
            
            for i in range(1, numberOfVols + 1):
                inputVol = locationToXmipp(i, self.inputFn)
                outputVol = locationToXmipp(i, self.outputStk)
                localArgs = self._segmentLocalArgs(self.inputFn, fnMask, args)
                maskArgs = self._segMentMaskArgs(inputVol, outputVol, fnMask)
                self._segmentVolume(localArgs, maskArgs, fnMask)
    
    def _segmentVolume(self, localArgs, maskArgs, fnMask):
        self.runJob("xmipp_volume_segment", localArgs)
        if exists(fnMask):
            self.runJob("xmipp_transform_mask", maskArgs)
            cleanPath(fnMask)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _argsChangeHand(self):
        if self.isFirstStep:
            if self._isSingleInput():
                args = "-i %s -o %s" % (self.inputFn, self.outputStk)
            else:
                args = "-i %s -o %s --save_metadata_stack %s --keep_input_columns" % (self.inputFn, self.outputStk, self.outputMd)
        else:
            args = "-i %s" % self.outputStk
        args += " --flipX"
        return args
    
    def _argsRandomize(self):
        if self.isFirstStep:
            if self._isSingleInput():
                args = "-i %s -o %s" % (self.inputFn, self.outputStk)
            else:
                args = "-i %s -o %s --save_metadata_stack %s --keep_input_columns" % (self.inputFn, self.outputStk, self.outputMd)
        else:
            args = "-i %s" % self.outputStk
        
        samplingRate = self.inputVolumes.get().getSamplingRate()
        resol = self.maxResolutionRandomize.get()
        args += " --freq continuous %f %f" % (float(resol),samplingRate)
        return args
    
    def _argsSymmetrize(self):
        if self.isFirstStep:
            if self._isSingleInput():
                args = "-i %s -o %s" % (self.inputFn, self.outputStk)
            else:
                args = "-i %s -o %s --save_metadata_stack %s --keep_input_columns" % (self.inputFn, self.outputStk, self.outputMd)
        else:
            args = "-i %s" % self.outputStk
        
        symmetry = self.symmetryGroup.get()
        symmetryAggregation = self.aggregation.get()
        
        # Validation done in the _validate function
#         if symmetry != 'c1':
        args += " --sym %s" % symmetry
        
        if symmetryAggregation == "sum":
            args += " --sum"
                
        return args
    
    def _argsAdjust(self):
        # ToDo Fix
        projFn = self.setOfProjections.get()
        args = " -m %s" % projFn
        return args
    
    def _adjustLocalArgs(self, inputVol, outputVol, args):
            if self.isFirstStep:
                localArgs = "-i %s -o %s" % (inputVol, outputVol) + args
            else:
                localArgs = "-i %s" % outputVol + args
            return localArgs
    
    def _argsSegment(self):
        segmentationType = self.segmentationType.get()
        segmentationMass = self.segmentationMass.get()
        ts = self._getSize()
        
        args = " --method "
        if segmentationType == "Voxel mass":
            args += "voxel_mass %d" % (int(segmentationMass))
        elif segmentationType == "Aminoacid mass":
            args += "aa_mass %d %f" % (int(segmentationMass),float(ts))
        elif segmentationType == "Dalton mass":
            args += "dalton_mass %d %f" % (int(segmentationMass),float(ts))
        else:
            args += "otsu"
        
        return args
    
    def _segmentLocalArgs(self, inputVol, fnMask, args):
        return "-i %s -o %s " % (inputVol, fnMask) + args
    
    def _segMentMaskArgs(self, inputVol, outputVol, fnMask):
            if self.isFirstStep:
                maskArgs = "-i %s -o %s" % (inputVol, outputVol)
            else:
                maskArgs = "-i %s" % outputVol
            maskArgs += " --mask binary_file %s" % fnMask
            return maskArgs
    
    def _argsNormalize(self):
        if self.isFirstStep:
            if self._isSingleInput():
                args = "-i %s -o %s" % (self.inputFn, self.outputStk)
            else:
                args = "-i %s -o %s --save_metadata_stack %s --keep_input_columns" % (self.inputFn, self.outputStk, self.outputMd)
        else:
            args = "-i %s" % self.outputStk
        
        maskRadius = self.backRadius.get()
        if maskRadius <= 0:
            size = self._getSize()
            maskRadius = size/2
        
        args += " --method NewXmipp --background circle %d" % int(maskRadius)
        return args
    
    def _argsInvert(self):
        if self.isFirstStep:
            if self._isSingleInput():
                args = "-i %s -o %s" % (self.inputFn, self.outputStk)
            else:
                args = "-i %s -o %s --save_metadata_stack %s --keep_input_columns" % (self.inputFn, self.outputStk, self.outputMd)
        else:
            args = "-i %s" % self.outputStk
        args += self._argsCommonInvert()
        return args
    
    def _argsThreshold(self):
        if self.isFirstStep:
            if self._isSingleInput():
                args = "-i %s -o %s" % (self.inputFn, self.outputStk)
            else:
                args = "-i %s -o %s --save_metadata_stack %s --keep_input_columns" % (self.inputFn, self.outputStk, self.outputMd)
        else:
            args = "-i %s" % self.outputStk
        args += self._argsCommonThreshold()
        return args

    def _validate(self):
        validateMsgs = []
        
        if not self.inputVolumes.hasValue():
            validateMsgs.append('Please provide an initial volume(s).')
            
        if self.doNormalize.get():
            size = int(self._getSize()/2)
            
            if self.backRadius.get() > size:
                validateMsgs.append('Set a valid Background radius less than %d' % size)
        
        if self.doSymmetrize.get():
            if self.symmetryGroup.get() == 'c1':
                validateMsgs.append('c1 is not a valid symmetry group.\n If you do not want to symmetrize set the field Symmetrize to not.')
                
        return validateMsgs
    
    def _summary(self):
        summary = []
        summary.append("Input volumes:  %s" % self.inputVolumes.get().getNameId())
        
        if not hasattr(self, 'outputVol'):
            summary.append("Output volumes not ready yet.")
        else:
            summary.append("Output volumes: %s" % self.outputVol)
        
        return summary
    
    def _methods(self):
        return self._summary()

    #--------------------------- UTILS functions ---------------------------------------------------
    def _getSize(self):
        """ get the size of either Volume or SetOfVolumes object"""
        if isinstance(self.inputVolumes.get(), Volume):
            Xdim = self.inputVolumes.get().getDim()[0]
        else:
            Xdim = self.inputVolumes.get().getDimensions()[0]
        return Xdim
