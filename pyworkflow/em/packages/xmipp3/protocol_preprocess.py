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
from protocol_process import XmippProcess, XmippProcessParticles, XmippProcessVolumes
from pyworkflow.utils.path import cleanPath
from pyworkflow.em.constants import *
from constants import *
from convert import locationToXmipp



class XmippPreprocess():
    """ This class has the common functions to preprocess objects like SetOfParticles, Volume or SetOfVolumes. """
    
    def __init__(self, **args):
        pass
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
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
    def _insertCommonSteps(self, outputFn):
        if self.doInvert:
            self._insertFunctionStep("invertStep", outputFn)
        
        if self.doThreshold:
            self._insertFunctionStep("thresholdStep", outputFn, self.threshold.get())
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def invertStep(self, outputFn):
        args = ' -i %s --mult -1' % outputFn
        self.runJob('xmipp_image_operate', args)
    
    def thresholdStep(self, outputFn, threshold):
        args = " -i %(outputFn)s --select below %(threshold)f "
        fillStr = self.getEnumText('fillType')
        args += "--substitute %(fillStr)s "
        
        if self.fillType == MASK_FILL_VALUE:
            args += " %f" % self.fillValue.get()
            
        self.runJob("xmipp_transform_threshold", args % locals())


class XmippProtPreprocessParticles(XmippProcessParticles, XmippPreprocess):
    """ Apply some filter to SetOfParticles """
    _label = 'preprocess particles'
    
    def __init__(self, **args):
        XmippProcessParticles.__init__(self, **args)
        XmippPreprocess.__init__(self, **args)
    
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
                      help='OldXmipp (mean(Image)=0, stddev(Image)=1).  \n  '
                           'NewXmipp (mean(background)=0, stddev(background)=1)  \n  '
                           'Ramp (subtract background+NewXmipp).  \n  ',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('backRadius', IntParam, default=-1, condition='doNormalize',
                      label='Background radius',
                      help='Pixels outside this circle are assumed to be noise and their stddev '
                      'is set to 1. Radius for background circle definition (in pix.). '
                      'If this value is 0, then half the box size is used.', 
                      expertLevel=LEVEL_ADVANCED)
        XmippPreprocess._defineProcessParams(self, form)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertProcessStep(self, inputFn, outputFn, outputMd):
        if self.doRemoveDust:
            self._insertFunctionStep("removeDustStep", outputFn)
        
        if self.doNormalize:
            self._insertFunctionStep("normalizeStep", outputFn, self.normType.get(), self.backRadius.get())
        
        self._insertCommonSteps(outputFn)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def removeDustStep(self, outputFn):
        threshold = self.thresholdDust.get()
        self.runJob('xmipp_transform_filter','-i %(outputFn)s --bad_pixels outliers %(threshold)f' % locals())
    
    def normalizeStep(self, ImagesMd, normType, bgRadius):
        args = "-i %(ImagesMd)s "
        radii = self._getSize()
        print "RADIUS:", radii
        if bgRadius <= 0:
            bgRadius = int(radii)
        
        if normType == "OldXmipp":
            args += "--method OldXmipp"
        elif normType == "NewXmipp":
            args += "--method NewXmipp --background circle %(bgRadius)d"
        else:
            args += "--method Ramp --background circle %(bgRadius)d"
        self.runJob("xmipp_transform_normalize", args % locals())
    
    def convertStep(self):
        """ convert if necessary"""
        imgSet = self.inputParticles.get()
        imgSet.writeStack(self.outputStk)
    
    def createOutputStep(self):
        inImgSet = self.inputParticles.get()
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(inImgSet)
        
        for i, img in enumerate(inImgSet):
            j = i + 1 
            img.setLocation(j, self.outputStk)
            imgSet.append(img)
        
        self._defineOutputs(outputParticles=imgSet)
        self._defineTransformRelation(inImgSet, imgSet)
    
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
        summary.append("Input particles:  %s" % self.inputParticles.get().getNameId())
        
        if not hasattr(self, 'outputVol'):
            summary.append("Output particles not ready yet.")
        else:
            summary.append("Output particles: %s" % self.outputParticles.getNameId())
        
        return summary
    
    def _methods(self):
        return self._summary()
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _getSize(self):
        """ get the size of SetOfParticles object"""
        Xdim, _, _, _ = self.inputParticles.get().getDimensions()
        size = int(Xdim/2)
        return size


class XmippProtPreprocessVolumes(XmippProcessVolumes, XmippPreprocess):
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
        XmippProcessVolumes.__init__(self, **args)
        XmippPreprocess.__init__(self, **args)
    
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
        form.addParam('symmetryGroup', TextParam, default='c1',
                      label="Symmetry group", condition='doSymmetrize',
                      help='See [[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry][Symmetry]]'
                      'for a description of the symmetry groups format, If no symmetry is present, give c1.')
        form.addParam('aggregation', EnumParam, choices=['Average', 'Sum'], 
                      display=EnumParam.DISPLAY_COMBO,
                      default=0, label='Aggregation mode', condition = 'doSymmetrize',
                      help='Symmetrized volumes can be averaged or summed.')
        # Adjust gray values
        form.addParam('doAdjust', BooleanParam, default=False,
                      label="Adjust gray values", 
                      help='Adjust input gray values so that it is compatible with a set of projections.')
        form.addParam('setOfProjections', PointerParam, pointerClass='SetOfImages',
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
        # ToDo: add wizard
        form.addParam('backRadius', FloatParam, default=-1,
                      label="Mask Radius", condition='doNormalize',
                      help='In pixels. Set to -1 for half of the size of the volume.')
        XmippPreprocess._defineProcessParams(self, form)
        
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertProcessStep(self, inputFn, outputFn, outputMd):
        if self.doChangeHand:
            self._insertFunctionStep("changeHandStep", outputFn)
        
        if self.doRandomize:
            self._insertFunctionStep("randomizeStep", outputFn, self.maxResolutionRandomize.get())
        
        if self.doSymmetrize:
            self._insertFunctionStep("symmetrizeStep", outputFn, self.symmetryGroup.get(), self.aggregation.get())
        
        if self.doAdjust:
            self._insertFunctionStep("adjustStep", outputFn, self.setOfProjections.get())

        if self.doSegment:
            self._insertFunctionStep("segmentStep", outputFn, self.segmentationType.get(), self.segmentationMass.get())
        
        if self.doNormalize:
            self._insertFunctionStep("normalizeStep", outputFn, self.backRadius.get())
        
        self._insertCommonSteps(outputFn)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def changeHandStep(self, outModel):
        self.runJob("xmipp_transform_mirror","-i %s --flipX" % outModel)
    
    def randomizeStep(self, outModel, maxResolution):
        samplingRate = self.inputVolumes.get().getSamplingRate()
        self.runJob("xmipp_transform_randomize_phases", "-i %s --freq continuous %f %f" % (outModel,float(maxResolution),samplingRate))
    
    def symmetrizeStep(self, outModel, symmetry, symmetryAggregation):
        if symmetry!='c1':
            args="-i %s --sym %s"%(outModel, symmetry)
            if symmetryAggregation=="sum":
                args+=" --sum"
            self.runJob("xmipp_transform_symmetrize", args)
    
    def adjustStep(self, outModel, setOfImages):
        self.runJob("xmipp_adjust_volume_grey_levels", "-i %s -m %s" % (outModel, setOfImages))
    
    def segmentStep(self, outModel, segmentationType, segmentationMass):
        
        if self.singleVolume:
            self._segmentVolume(outModel, segmentationType, segmentationMass)
        else:
            numberOfParticles = self.inputVolumes.get().getSize()
            for i in range(1, numberOfParticles + 1):
                
                volOutModel = locationToXmipp(i, outModel)
                self._segmentVolume(volOutModel, segmentationType, segmentationMass)
        
    def _segmentVolume(self, outModel, segmentationType, segmentationMass):
        fnMask = self._getTmpPath("mask.vol")
        ts = self._getSize()
        args = "-i %s -o %s --method " % (outModel, fnMask)
        
        if segmentationType == "Voxel mass":
            args += "voxel_mass %d" % (int(segmentationMass))
        elif segmentationType == "Aminoacid mass":
            args += "aa_mass %d %f" % (int(segmentationMass),float(ts))
        elif segmentationType == "Dalton mass":
            args += "dalton_mass %d %f" % (int(segmentationMass),float(ts))
        else:
            args += "otsu"
        self.runJob("xmipp_volume_segment", args)
        
        if exists(fnMask):
            self.runJob("xmipp_transform_mask", "-i %s --mask binary_file %s" % (outModel, fnMask))
            cleanPath(fnMask)
    
    def normalizeStep(self, outModel, maskRadius):
        if maskRadius <= 0:
            size = self._getSize()
            maskRadius = size/2
        self.runJob("xmipp_transform_normalize", "-i %s --method NewXmipp --background circle %d" % (outModel, int(maskRadius)))
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        validateMsgs = []
        
        if not self.inputVolumes.hasValue():
            validateMsgs.append('Please provide an initial volume(s).')
            
        if self.doNormalize.get():
            size = int(self._getSize()/2)
            
            if self.backRadius.get() > size:
                validateMsgs.append('Set a valid Background radius less than %d' % size)
            
        return validateMsgs
    
    def _summary(self):
        summary = []
        summary.append("Input volumes:  %s" % self.inputVolumes.get().getNameId())
        
        if not hasattr(self, 'outputVol'):
            summary.append("Output volumes not ready yet.")
        else:
            summary.append("Output volumes: %s" % self.outputVol)
        
        return summary
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _getSize(self):
        """ get the size of either Volume or SetOfVolumes object"""
        if isinstance(self.inputVolumes.get(), Volume):
            Xdim, _, _, _ = self.inputVolumes.get().getDim()
        else:
            Xdim, _, _, _ = self.inputVolumes.get().getDimensions()
        return Xdim
