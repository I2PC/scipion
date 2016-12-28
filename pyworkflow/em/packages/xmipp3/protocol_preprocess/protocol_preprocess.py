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

from pyworkflow.em import *
from pyworkflow.utils import *  
from pyworkflow.protocol.params import *
from protocol_process import XmippProcessParticles, XmippProcessVolumes
from pyworkflow.utils.path import cleanPath
from ..constants import *
from pyworkflow.em.packages.xmipp3.convert import getImageLocation
from ..convert import locationToXmipp, writeSetOfParticles


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
            protocol._insertFunctionStep("invertStep", args, changeInserts)

        if protocol.doThreshold:
            args = protocol._argsThreshold()
            protocol._insertFunctionStep("thresholdStep", args, changeInserts)

    #--------------------------- UTILS functions ---------------------------------------------------
    @classmethod
    def _argsCommonInvert(cls):
        args = ' --mult -1'
        return args
    
    @classmethod
    def _argsCommonThreshold(cls, protocol):
        args = " --select %s %f" % (protocol.getEnumText('thresholdType'), protocol.threshold)
        fillStr = protocol.getEnumText('fillType')
        args += " --substitute %s " % fillStr
        
        if protocol.fillType == MASK_FILL_VALUE:
            args += " %f" % protocol.fillValue
        return args
    

class XmippProtPreprocessParticles(XmippProcessParticles):
    """ Preprocess a set of particles. You can remove dust, normalize, apply threshold, etc """
    _label = 'preprocess particles'

    # Automatic Particle rejection enum
    REJ_NONE = 0
    REJ_MAXZSCORE = 1
    REJ_PERCENTAGE =2
    
    def __init__(self, **kwargs):
        XmippProcessParticles.__init__(self, **kwargs)
    
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
        form.addParam('doCenter', BooleanParam, default=False,
                      label='Center images')
        XmippPreprocessHelper._defineProcessParams(form)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertProcessStep(self):
        self.isFirstStep = True
        # this is for when the options selected has changed and the protocol is resumed
        changeInserts = [self.doRemoveDust, self.doNormalize, self.doInvert, self.doThreshold, self.doCenter]
        
        if self.doRemoveDust:
            args = self._argsRemoveDust()
            self._insertFunctionStep("removeDustStep", args, changeInserts)
        
        if self.doNormalize:
            args = self._argsNormalize()
            self._insertFunctionStep("normalizeStep", args, changeInserts)
        
        if self.doCenter:
            args = self._argsCenter()
            self._insertFunctionStep("centerStep", args, changeInserts)
        
        XmippPreprocessHelper._insertCommonSteps(self, changeInserts)
        
    #--------------------------- STEPS functions ---------------------------------------------------
    def invertStep(self, args, changeInserts):
        self.runJob('xmipp_image_operate', args)

    def thresholdStep(self, args, changeInserts):
        self.runJob("xmipp_transform_threshold", args)
        
    def removeDustStep(self, args, changeInserts):
        self.runJob('xmipp_transform_filter', args)
    
    def normalizeStep(self, args, changeInserts):
        self.runJob("xmipp_transform_normalize", args)
    
    def centerStep(self, args, changeInserts):
        self.runJob("xmipp_transform_center_image", args % locals())
    
    def sortImages(self, outputFn, outputMd):
        pass

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
        if hasattr(self, 'outputParticles'):
            methods.append("Input particles %s of %s elements" % (self.getObjectTag('inputParticles'), self.inputParticles.get().getSize()))
            if self.doNormalize:
                methods.append("The background was normalized with %s method." % self.getEnumText('normType'))
            if self.doInvert:
                methods.append("The contrast was inverted")
            if self.doThreshold:
                methods.append("Pixels with values below %f was removed" % self.threshold.get())
            methods.append('Output set: %s'%self.getObjectTag('outputParticles'))
        return methods
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _argsRemoveDust(self):
        if self.isFirstStep:
            args = "-i %s -o %s --save_metadata_stack %s --keep_input_columns" % (self.inputFn, self.outputStk, self.outputMd)
            self._setFalseFirstStep()
        else:
            args = "-i %s" % self.outputStk
        args += " --bad_pixels outliers %f" % self.thresholdDust.get()
        return args
    
    def _argsNormalize(self):
        if self.isFirstStep:
            args = "-i %s -o %s --save_metadata_stack %s --keep_input_columns" % (self.inputFn, self.outputStk, self.outputMd)
            self._setFalseFirstStep()
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
            self._setFalseFirstStep()
        else:
            args = "-i %s" % self.outputStk
        args += XmippPreprocessHelper._argsCommonInvert()
        return args
    
    def _argsThreshold(self):
        if self.isFirstStep:
            args = "-i %s -o %s --save_metadata_stack %s --keep_input_columns" % (self.inputFn, self.outputStk, self.outputMd)
            self._setFalseFirstStep()
        else:
            args = "-i %s" % self.outputStk
        args += XmippPreprocessHelper._argsCommonThreshold(self)
        return args
    
    def _argsCenter(self):
        if self.isFirstStep:
            args = "-i %s -o %s --save_metadata_stack %s" % (self.inputFn, self.outputStk, self.outputMd)
            self._setFalseFirstStep()
        else:
            args = "-i %s" % self.outputStk
        return args
    
    def _getSize(self):
        """ get the size of SetOfParticles object"""
        Xdim = self.inputParticles.get().getDimensions()[0]
        size = int(Xdim/2)
        return size
    
    def _setFalseFirstStep(self):
        if self.isFirstStep:
                self.isFirstStep = False


class XmippProtPreprocessVolumes(XmippProcessVolumes):
    """ Protocol for Xmipp-based preprocess for volumes """
    import pyworkflow.em.metadata as md
    
    _label = 'preprocess volumes'
    
    # Aggregation constants
    AGG_AVERAGE=0
    AGG_SUM=1
    
    # Segmentation type
    SEG_VOXEL=0
    SEG_AMIN=1
    SEG_DALTON=2
    SEG_AUTO=3


    def __init__(self, **kwargs):
        XmippProcessVolumes.__init__(self, **kwargs)
    
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
        form.addParam('volumeMask', PointerParam, pointerClass='VolumeMask', allowsNull=True,
                      label='Mask volume', condition='doSymmetrize'
                      )
        form.addParam('doWrap', BooleanParam, default=True,
                      label="Wrap", condition='doSymmetrize',
                      help='by default, the image/volume is wrapped')
        # Adjust gray values
        form.addParam('doAdjust', BooleanParam, default=False,
                      label="Adjust gray values", 
                      help='Adjust input gray values so that it is compatible with a set of projections.') 
        form.addParam('inputImages', PointerParam, pointerClass='SetOfParticles, SetOfAverages, SetOfClasses2D',
                      label="Set of particles", condition='doAdjust',
                      help='Set of images to which the model should conform. The set of images should have the'
                      'final pixel size and the final size of the model.')
        form.addParam('sigSymGroup', TextParam, default='c1',
                      label="Symmetry group", condition='doAdjust',
                      help='See [[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry][Symmetry]]'
                      'for a description of the symmetry groups format, If no symmetry is present, give c1.')  
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
        XmippPreprocessHelper._defineProcessParams(form)

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertProcessStep(self):
        self.isFirstStep = True
        # this is for when the options selected has changed and the protocol is resumed
        changeInserts = [self.doChangeHand, self.doRandomize, self.doSymmetrize, 
                         self.doAdjust, self.doSegment, self.doInvert, self.doNormalize, 
                         self.doThreshold]

        if self.doChangeHand:
            args = self._argsChangeHand()
            self._insertFunctionStep("changeHandStep", args, changeInserts)
        
        if self.doRandomize:
            args = self._argsRandomize()
            self._insertFunctionStep("randomizeStep", args, changeInserts)
        
        if self.doSymmetrize:
            args = self._argsSymmetrize()
            self._insertFunctionStep("symmetrizeStep", args, changeInserts)
        
        if self.doAdjust:
            self._insertFunctionStep("projectionStep", changeInserts)
            self._insertFunctionStep("adjustStep", self.isFirstStep, changeInserts)
            if self.isFirstStep:
                self.isFirstStep = False

        if self.doSegment:
            args = self._argsSegment()
            self._insertFunctionStep("segmentStep", args, changeInserts)
            if self.isFirstStep:
                self.isFirstStep = False
        
        if self.doNormalize:
            args = self._argsNormalize()
            self._insertFunctionStep("normalizeStep", args, changeInserts)
        
        XmippPreprocessHelper._insertCommonSteps(self, changeInserts)

    #--------------------------- STEPS functions ---------------------------------------------------
    def invertStep(self, args, changeInserts):
        self.runJob('xmipp_image_operate', args)
    
    def thresholdStep(self, args, changeInserts):
        self.runJob("xmipp_transform_threshold", args)
        
    def removeDustStep(self, args, changeInserts):
        self.runJob('xmipp_transform_filter', args)
    
    def normalizeStep(self, args, changeInserts):
        self.runJob("xmipp_transform_normalize", args)
        
    def changeHandStep(self, args, changeInserts):
        self.runJob("xmipp_transform_mirror", args)
    
    def randomizeStep(self, args, changeInserts):
        self.runJob("xmipp_transform_randomize_phases", args)
    
    def symmetrizeStep(self, args, changeInserts):
        self.runJob("xmipp_transform_symmetrize", args)
    
    def projectionStep(self, changeInserts):
        partSet = self.inputImages.get()
        imgsFn = self._getTmpPath('input_images.xmd')
        
        if partSet.getSize() > 200:
            newPartSet = self._getRandomSubset(partSet, 200)
        else:
            newPartSet = partSet
            
        writeSetOfParticles(newPartSet, imgsFn)
        
        if not partSet.hasAlignmentProj():
            params = {'imgsFn': imgsFn,
                      'dir': self._getTmpPath(),
                      'vols': self.inputFn,
                      'symmetryGroup': self.sigSymGroup.get(),
                      }
            sigArgs = '-i %(imgsFn)s --initvolumes %(vols)s --odir %(dir)s --sym %(symmetryGroup)s'\
            ' --alpha0 0.005 --dontReconstruct' % params
            self.runJob("xmipp_reconstruct_significant", sigArgs)
    
    def adjustStep(self, isFirstStep, changeInserts):
        if isFirstStep:
            inputFn = self.inputFn
        else:
            if self._isSingleInput():
                inputFn = self.outputStk
            else:
                inputFn = self.outputMd

        if self._isSingleInput():
            args = self._argsAdjust(0)
            localArgs = self._adjustLocalArgs(inputFn, self.outputStk, args)
            self.runJob("xmipp_transform_adjust_volume_grey_levels", localArgs)
        else:
            volMd = md.MetaData(self.inputFn)
            outVolMd = md.MetaData(self.outputMd)
            for objId in volMd:
                args = self._argsAdjust(objId-1)
                inputVol = volMd.getValue(md.MDL_IMAGE, objId)
                outputVol = outVolMd.getValue(md.MDL_IMAGE, objId)
                localArgs = self._adjustLocalArgs(inputVol, outputVol, args)
                self.runJob("xmipp_transform_adjust_volume_grey_levels", localArgs)
    
    def segmentStep(self, args, changeInserts):
        fnMask = self._getTmpPath("mask.vol")
        if self.isFirstStep:
            inputFn = self.inputFn
        else:
            if self._isSingleInput():
                inputFn = self.outputStk
            else:
                inputFn = self.outputMd
        
        if self._isSingleInput():
            localArgs = self._segmentLocalArgs(inputFn, fnMask, args)
            maskArgs = self._segMentMaskArgs(inputFn, self.outputStk, fnMask)
            self._segmentVolume(localArgs, maskArgs, fnMask)
        else:
            volMd = md.MetaData(inputFn)
            outVolMd = md.MetaData(self.outputMd)
            for objId in volMd:
                inputVol = volMd.getValue(md.MDL_IMAGE, objId)
                outputVol = outVolMd.getValue(md.MDL_IMAGE, objId)
                localArgs = self._segmentLocalArgs(inputVol, fnMask, args)
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
            self._setFalseFirstStep()
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
            self._setFalseFirstStep()
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
            self._setFalseFirstStep()
        else:
            args = "-i %s" % self.outputStk

        symmetry   = self.symmetryGroup.get()
        doWrap     = self.doWrap.get()
        if self.volumeMask.get() is not None:
            fnVolumeMask = self.volumeMask.get().getFileName()
            doVolumeMask = True
        else:
            doVolumeMask = False

        ###########FILEFILEFILE
        symmetryAggregation = self.aggregation.get()

        # Validation done in the _validate function
        #         if symmetry != 'c1':
        args += " --sym %s " % symmetry

        if symmetryAggregation == "sum":
            args += " --sum"

        if not doWrap:
            args += " --dont_wrap "

        if doVolumeMask:
            if exists(fnVolumeMask):
                args += " --mask_in %s "%fnVolumeMask
            else:
                print('Error: mask %s does not exists'%fnVolumeMask)
        return args
    
    def _argsAdjust(self, number):
        if self.inputImages.get().hasAlignmentProj():
            fn = "input_images.xmd"
        else:
            fn = "images_iter001_%02d.xmd" % number
        args = " -m %s" % self._getTmpPath(fn)
        return args
    
    def _adjustLocalArgs(self, inputVol, outputVol, args):
            localArgs = "-i %s -o %s" % (inputVol, outputVol) + args
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
        print "self.isFirstStep, ", self.isFirstStep
        if self.isFirstStep:
            maskArgs = "-i %s -o %s" % (inputVol, outputVol)
            self._setFalseFirstStep()
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
            self._setFalseFirstStep()
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
            self._setFalseFirstStep()
        else:
            args = "-i %s" % self.outputStk
        args += XmippPreprocessHelper._argsCommonInvert()
        return args
    
    def _argsThreshold(self):
        if self.isFirstStep:
            if self._isSingleInput():
                args = "-i %s -o %s" % (self.inputFn, self.outputStk)
            else:
                args = "-i %s -o %s --save_metadata_stack %s --keep_input_columns" % (self.inputFn, self.outputStk, self.outputMd)
        else:
            args = "-i %s" % self.outputStk
        args += XmippPreprocessHelper._argsCommonThreshold(self)
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
    
    def _setFalseFirstStep(self):
        if self.isFirstStep:
                self.isFirstStep = False

    def _getRandomSubset(self, imgSet, numOfParts):
        if isinstance(imgSet, SetOfClasses2D):
            partSet = self._createSetOfParticles("_averages")
            for i, cls in enumerate(imgSet):
                img = cls.getRepresentative()
                img.setSamplingRate(cls.getSamplingRate())
                img.setObjId(i+1)
                partSet.append(img)
        else:
            partSet = imgSet
        
        if partSet.getSize() > numOfParts:
            newPartSet = SetOfParticles(filename=self._getTmpPath("particles_tmp.sqlite"))
            counter = 0
            for part in partSet.iterItems(orderBy='RANDOM()', direction='ASC'):
                if counter < numOfParts:
                    newPartSet.append(part)
                    counter =+ 1
                else:
                    break
        else:
            newPartSet = partSet
        
        return newPartSet
    
