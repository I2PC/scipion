# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
from math import floor
from xmipp3 import XmippProtocol
from convert import createXmippInputVolumes, readSetOfVolumes

#from xmipp3 import XmippProtocol

class XmippProtPreprocessVolumes(ProtPreprocessVolumes, XmippProtocol):
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


    def _defineParams(self, form):
        form.addSection(label='Input')
        # Volumes to process
        form.addParam('inputVolumes', PointerParam, label="Input Volume(s):", important=True, 
                      pointerClass='SetOfVolumes, Volume',
                      help='Can be a density volume or a SetOfVolumes.')
#         form.addParam('doResize', BooleanParam, default=False,
#                       label="Resize the volume(s)", 
#                       help='If set to True, resize the volume(s).')
#         form.addParam('outputVoxel', FloatParam, default=0, condition='doResize', label='Final voxel size (A/voxel):',
#                       help='Set a desire sampling rate.\n '
#                       'If Final box size is larger than 0, the sampling rate will be\n'
#                       'calculated by automatic estimation. If both Final voxel size and\n '
#                       'Final box size has valid values, the protocol ignores the Final box size\n'
#                       'than Final box size')
#         form.addParam('size', IntParam, default=0, condition='doResize',
#                       label='Final box size (voxels):', help='Set to 0 (or negative) for automatic estimation\n'
#                       'if a valid value of final voxel size is set. If both Final voxel size and\n'
#                       'Final box size has valid values, the protocol ignores the Final box size\n'
#                       'than Final box size')
        
        form.addSection(label='Preprocess steps')
        form.addParam('doWindow', BooleanParam, default=False,
                      label="Pad or crop volume(s)",
                      help='If set to True, resize the volume(s).')
        form.addParam('finalSize', IntParam, default=256, condition='doWindow',
                      label='Final dimensions (voxels):',
                      help='Set the desire value to modify the dimensions.\n'
                      'Important: With this transformation you NOT modify the sampling rate\n'
                      'of your SetOfVolumes (or volume). Only modify the dimensions')
        # Change hand
        form.addParam('doChangeHand', BooleanParam, default=False,
                      label="Change hand", 
                      help='Change hand by applying a mirror along X.')
        
        form.addParam('doRandomize', BooleanParam, default=False,
                      label="Randomize phases", 
                      help='Randomize phases beyond a certain frequency.')
        # Maximum Resolution
        form.addParam('maxResolutionRandomize', FloatParam, default=40,
                      label="Maximum Resolution", condition='doRandomize',
                      help='Angstroms.')
        # Maximum Resolution
        form.addParam('doFilter', BooleanParam, default=False,
                      label="Low pass filter", 
                      help='Low pass filter the volume.')
        form.addParam('maxResFilt', IntParam, default=40,
                      label="Maximum Resolution", condition='doFilter',
                      help='Angstroms.')        
        # Symmetry group
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
        # Spherical Mask
        form.addParam('doMask', BooleanParam, default=False,
                      label="Spherical Mask", 
                      help='Remove all pixels beyond a certain radius.')
        form.addParam('maskRadius', FloatParam, default=-1,
                      label="Maximum Resolution", condition='doMask',
                      help='In pixels. Set to -1 for half of the size of the volume.') 
        # Adjust gray values
        form.addParam('doAdjust', BooleanParam, default=False,
                      label="Adjust gray values", 
                      help='Adjust input gray values so that it is compatible with a set of projections.')
        form.addParam('setOfProjections', PointerParam, pointerClass='SetOfImages',
                      label="Set of projections", condition='doAdjust',
                      help='Metadata with a set of images to which the model should conform. The set of images should have the'
                      'final pixel size and the final size of the model.')
        # Normalize background
        form.addParam('doNormalize', BooleanParam, default=False,
                      label="Normalize background", 
                      help='Set background to have zero mean and standard deviation 1.')
        form.addParam('maskRadiusNormalize', FloatParam, default=-1,
                      label="Mask Radius", condition='doNormalize',
                      help='In pixels. Set to -1 for half of the size of the volume.')
        # Threshold
        form.addParam('doThreshold', BooleanParam, default=False,
                      label="Threshold", 
                      help='Remove voxels below a certain value.')
        form.addParam('threshold', FloatParam, default=0,
                      label="Threshold value", condition='doThreshold',
                      help='Grey value below which all voxels should be set to 0.')
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
        
#        form.addParallelSection(threads=1, mpi=8)         
    
    def _insertAllSteps(self):
        """The preprocess is realized for a set of volumes or a single volume"""
        
        # Convert SetOfVolumes to Xmipp Metadata
        volSet = self.inputVolumes.get()
        self.outModel = self._getPath("volumes.stk")
        samplingRate = volSet.getSamplingRate()
        finalSamplingRate = samplingRate

#         xDim, _, _, self.ndim = volSet.getDimensions()
#         if self.doResize.get():
# 
#             if self.outputVoxel.get() >> 0:
#                 finalSamplingRate = self.outputVoxel.get()
#                 factor = floor(samplingRate / finalSamplingRate)
#                 args = "-i %s -o %s --factor %f" % (self.inModel, self.outModel, factor)
#             elif self.size.get() >= 1:
#                 dimensions = self.size.get()
#                 finalSamplingRate = floor(self.size.get() * samplingRate / xDim)
#                 args = "-i %s -o %s --dim %d" % (self.inModel, self.outModel, dimensions)
#             
#             if not singleVolume:
#                 args+=" --save_metadata_stack"
#             self._insertRunJobStep(self, "xmipp_image_resize", args)
        
        self._insertFunctionStep('copyFilesStep', self.outModel)
        if self.doWindow:
            self._insertFunctionStep("windowStep", self.outModel, self.finalSize.get())
        
        if self.doChangeHand:
            self._insertFunctionStep("changeHandStep", self.outModel)
        
        if self.doRandomize:
            self._insertFunctionStep("randomizeStep", self.outModel, finalSamplingRate, self.maxResolutionRandomize.get())
        
        if self.doFilter:
            self._insertFunctionStep("filterStep", self.outModel, self.finalSize.get(), self.maxResFilt.get())
        
        if self.doSymmetrize:
            self._insertFunctionStep("symmetrizeStep", self.outModel, self.symmetryGroup.get(), self.aggregation.get())
        
        if self.doMask:
            self._insertFunctionStep("maskStep", self.outModel, self.maskRadius.get())
        
        if self.doAdjust:
            self._insertFunctionStep("adjust", self.outModel, self.setOfProjections.get())

        if self.doNormalize:
            self._insertFunctionStep("normalize", self.outModel, self.maskRadiusNormalize.get())
        
        if self.doThreshold:
            self._insertFunctionStep("thresholdFunc", self.outModel, self.threshold.get())
        
        if self.doSegment:
            self._insertFunctionStep("segment", self.outModel, self.segmentationType.get(), self.segmentationMass.get(),
                            self.finalSize.get())
        self._insertFunctionStep('createOutput', finalSamplingRate)
    
    def copyFilesStep(self, outModel):
    	""" Prepare the files to process """
    	volSet = self.inputVolumes.get()
		# Check volsMd is a volume or a stack
    	if isinstance(volSet, Volume):
			self.inModel  = volSet.getFileName()
			self.singleVolume = True
			args = "-i %s -o %s -t vol" % (self.inModel, outModel)
			self.runJob('xmipp_image_convert', args)
    	else:
    		volSet.writeStack(outModel)
    		self.singleVolume = False
    
    def windowStep(self, outModel, size):
        self.runJob("xmipp_transform_window", "-i %s --size %d" % (outModel, int(size)))
    
    def randomizeStep(self, outModel, ts, maxResolution):
        self.runJob("xmipp_transform_randomize_phases","-i %s --freq continuous %f %f"%(outModel,float(maxResolution),float(ts)))
    
    def filterStep(self, outModel, ts, maxResolution):
        self.runJob("xmipp_transform_filter","-i %s --fourier low_pass %f --sampling %f"%(outModel,float(maxResolution),float(ts)))
    
    def maskStep(self, outModel, maskRadius):
        if maskRadius==-1:
            md= xmipp.MetaData(outModel)
            Xdim, _, _, _, _ = xmipp.MetaDataInfo(md)    
            maskRadius=Xdim/2
        self.runJob("xmipp_transform_mask", "-i %s --mask circular %d"%(outModel,-int(maskRadius)))
    
    def symmetrizeStep(self, outModel, symmetry, symmetryAggregation):
        if symmetry!='c1':
            args="-i %s --sym %s"%(outModel, symmetry)
            if symmetryAggregation=="sum":
                args+=" --sum"
            self.runJob("xmipp_transform_symmetrize", args)
    
    def adjust(self, outModel, setOfImages):
        self.runJob("xmipp_adjust_volume_grey_levels", "-i %s -m %s"%(outModel, setOfImages))
    
    def normalize(self, outModel, maskRadius):
        if maskRadius==-1:
            md= xmipp.MetaData(outModel)
            Xdim, _, _, _, _ = xmipp.MetaDataInfo(md)    
            maskRadius=Xdim/2
        self.runJob("xmipp_transform_normalize", "-i %s --method NewXmipp --background circle %d" % (outModel, int(maskRadius)))
    
    def thresholdFunc(self, outModel, threshold):
        self.runJob("xmipp_transform_threshold", "-i %s --select below %f --substitute value 0" % (outModel, float(threshold)))
    
    def changeHandStep(self, outModel):
        self.runJob("xmipp_transform_mirror","-i %s --flipX"%outModel)
    
    def segment(self, outModel, segmentationType, segmentationMass, ts):
        fnMask=self._getPath("mask.vol")
        args="-i %s -o %s --method "%(outModel, fnMask)
        
        if segmentationType == "Voxel mass":
            args += "voxel_mass %d" % (int(segmentationMass))
        elif segmentationType == "Aminoacid mass":
            args += "aa_mass %d %f" % (int(segmentationMass),float(ts))
        elif segmentationType == "Dalton mass":
            args += "dalton_mass %d %f" % (int(segmentationMass),float(ts))
        else:
            args += "otsu"
        
        self.runJob("xmipp_volume_segment",args)
        
        if exists(fnMask):
            self.runJob("xmipp_transform_mask", "-i %s --mask binary_file %s" % (outModel, fnMask))
            
    def createOutput(self, samplingRate):
        
        if self.singleVolume:
            vol = Volume()
            vol.setFileName(self.outModel)
            self._defineOutputs(outputVol=vol)
        else:
            volumes = self._createSetOfVolumes()
            volumes.setSamplingRate(samplingRate)
            
            _, _, _, nDim = self.inputVolumes.get().getDimensions()
            for i in range(1,nDim+1):
                vol=Volume()
                vol.setLocation(i,self.outModel)
                volumes.append(vol)
            self._defineOutputs(outputVol=volumes)

        self._defineTransformRelation(self.inputVolumes.get(), self.outputVol)
    
    def _validate(self):
        validateMsgs = []
        
        if not self.inputVolumes.hasValue():
            validateMsgs.append('Please provide an initial volume(s).')
#         if self.doResize.get():
#             if self.outputVoxel.get() <= 0 and self.size.get() <= 0:
#             	validateMsgs.append('Please enter a valid value for Final voxel size or Final box size.')
        
        return validateMsgs
    
    def _summary(self):
        summary = []
        summary.append("Input volumes:  %s" % self.inputVolumes.get().getNameId())
        
        if not hasattr(self, 'outputVol'):
            summary.append("Output volumes not ready yet.")
        else:
            summary.append("Output volumes: %s" % self.outputVol)
        
        return summary
        
    