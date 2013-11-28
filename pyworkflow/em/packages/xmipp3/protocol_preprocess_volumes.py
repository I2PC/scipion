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
    _label = 'Xmipp Preprocess Volumes'
      
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
        form.addParam('inputVolumes', PointerParam, label="Volumes to process", important=True, 
                      pointerClass='SetOfVolumes',
                      help='Can be a density volume or a PDB file.')  
        form.addParam('inputVoxel', FloatParam, label='Input voxel size (A/voxel)')
        
        form.addSection(label='Preprocess steps')
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
                      help='See [http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry]'
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
        
        form.addSection(label='Output')
        form.addParam('outputVoxel', FloatParam, default=1.0, label='Final voxel size (A/voxel):')
        form.addParam('finalSize', FloatParam, default=1, label='Final box size (voxels):', help='Set to -1 for automatic estimation.')

#        form.addParallelSection(threads=1, mpi=8)         
    
    def _defineSteps(self):
        """The preprocess is realized for a set of volumes"""
        
        # Convert SetOfVolumes to Xmipp Metadata
        self.volsFn = createXmippInputVolumes(self, self.inputVolumes.get())
        self.volsMd = xmipp.MetaData(self.volsFn)
        
        # Check volsMd is a volume or a stack
        if self.volsMd.size()==1:
            self.outModel=self._getPath("volume.vol")
            self.singleVolume=True
        else:
            self.outModel=self._getPath("volumes.stk")
            self.singleVolume=False
        
        if self.inputVoxel.get() != self.outputVoxel.get() or self.finalSize.get()!=-1:
            if self.inputVoxel.get() != self.outputVoxel.get() and self.finalSize.get()==-1:
                x, _, _, _, _ = self.inputVolumes.get().getDimensions()
                fnSize = floor(x/(self.outputVoxel.get()/self.inputVoxel.get()))
                self.finalSize.set(fnSize)
            
            self._insertFunctionStep("changeSamplingRateAndOrBox",self.volsMd, self.outModel,
                            self.singleVolume, self.inputVoxel.get(),self.outputVoxel.get(),
                            self.finalSize.get())
        else:
            if self.singleVolume:
                copyFile(self.volsFn, self.outModel)
            else:
                self._insertRunJobStep('xmipp_image_convert', "-i %s -o %s --save_metadata_stack --track_origin"%(self.volsFn,self.outModel))
         
        if self.doChangeHand:
            self._insertFunctionStep("changeHand", self.outModel)
        
        if self.doRandomize:
            self._insertFunctionStep("randomize", self.outModel, self.finalSize.get(), self.maxResolutionRandomize)
        
        if self.doFilter:
            self._insertFunctionStep("filter", self.outModel, self.finalSize.get(), self.maxResFilt.get())
        
        if self.doSymmetrize:
            self._insertFunctionStep("symmetrize", self.outModel, self.symmetryGroup.get(), self.aggregation.get())
        
        if self.doMask:
            self._insertFunctionStep("mask", self.outModel, self.maskRadius.get())
        
        if self.doAdjust:
            self._insertFunctionStep("adjust", self.outModel, self.setOfProjections.get())

        if self.doNormalize:
            self._insertFunctionStep("normalize", self.outModel, self.MaskRadiusNormalize.get())
        
        if self.doThreshold:
            self._insertFunctionStep("threshold", self.outModel, self.threshold.get())
        
        if self.doSegment:
            self._insertFunctionStep("segment", self.outModel, self.segmentationType.get(), self.segmentationMass.get(),
                            self.finalSize.get())
         
         
    def window(self,input,outModel,singleVolume,size):
        if size>0:
            args="-i %s --size %d"%(input,int(size))
            if not singleVolume:
                args+=" --save_metadata_stack"
            if input!=outModel:
                args+=" -o %s"%outModel
            
            self._insertRunJobStep(self,"xmipp_transform_window",args)

    def scale(self,input,outModel,singleVolume,scale):
        args="-i %s -o %s --scale %f --dont_wrap"%(input,outModel,scale)
        if not singleVolume:
            args+=" --save_metadata_stack"
        
        self._insertRunJobStep(self,"xmipp_transform_geometry",args)
    
    def changeSamplingRateAndOrBox(self,inModel,outModel,singleVolume,initialTs,finalTs,size):
        if initialTs==finalTs:
            self.window(inModel,outModel,singleVolume,size)
        elif initialTs<finalTs:
            self.scale(inModel,outModel,singleVolume,initialTs/finalTs)
            self.window(outModel, outModel,singleVolume,size)
        else:
            self.window(inModel,outModel,singleVolume,size)
            self.scale(outModel,outModel,singleVolume,initialTs/finalTs)

    def randomize(self, outModel, ts, maxResolution):
        self._insertRunJobStep("xmipp_transform_randomize_phases","-i %s --freq continuous %f %f"%(outModel,float(maxResolution),float(ts)))
    
    def filter(self, outModel, ts, maxResolution):
        self._insertRunJobStep("xmipp_transform_filter","-i %s --fourier low_pass %f --sampling %f"%(outModel,float(maxResolution),float(ts)))
    
    def mask(self, outModel, maskRadius):
        if maskRadius==-1:
            md= xmipp.MetaData(outModel)
            Xdim, _, _, _, _ = xmipp.MetaDataInfo(md)    
            maskRadius=Xdim/2
        self._insertRunJobStep("xmipp_transform_mask","-i %s --mask circular %d"%(outModel,-int(maskRadius)))
    
    def symmetrize(self, outModel, symmetry, symmetryAggregation):
        if symmetry!='c1':
            args="-i %s --sym %s"%(outModel, symmetry)
            if symmetryAggregation=="sum":
                args+=" --sum"
            self._insertRunJobStep("xmipp_transform_symmetrize",args)
    
    def adjust(self, outModel, setOfImages):
        self._insertRunJobStep("xmipp_adjust_volume_grey_levels","-i %s -m %s"%(outModel, setOfImages))
    
    def normalize(self, outModel, maskRadius):
        if maskRadius==-1:
            md= xmipp.MetaData(outModel)
            Xdim, _, _, _, _ = xmipp.MetaDataInfo(md)    
            maskRadius=Xdim/2
        self._insertRunJobStep("xmipp_transform_normalize","-i %s --method NewXmipp --background circle %d"%(outModel, int(maskRadius)))
    
    def threshold(self, outModel, threshold):
        self._insertRunJobStep("xmipp_transform_threshold","-i %s --select below %f --substitute value 0"%(outModel, float(threshold)))
    
    def changeHand(self, outModel):
        self._insertRunJobStep("xmipp_transform_mirror","-i %s --flipX"%outModel)
    
    def segment(self, outModel, segmentationType, segmentationMass, ts):
        fnMask=self._getPath("mask.vol")
        args="-i %s -o %s --method "%(outModel, fnMask)
        
        if segmentationType== "Voxel mass":
            args+="voxel_mass %d"%(int(segmentationMass))
        elif segmentationType=="Aminoacid mass":
            args+="aa_mass %d %f"%(int(segmentationMass),float(ts))
        elif segmentationType=="Dalton mass":
            args+="dalton_mass %d %f"%(int(segmentationMass),float(ts))
        else:
            args+="otsu"
        
        self._insertRunJobStep("xmipp_volume_segment",args)
        
        if exists(fnMask):
            self._insertRunJobStep("xmipp_transform_mask","-i %s --mask binary_file %s"%(outModel, fnMask))
            
    def createOutput(self):
        pass
        
    def _validate(self):
        return []
    
    def _summary(self):
        pass

    