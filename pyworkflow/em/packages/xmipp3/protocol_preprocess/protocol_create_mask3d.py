# **************************************************************************
# *
# * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
#                J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from ..convert import getImageLocation
from ..constants import *
from geometrical_mask import *

SOURCE_VOLUME=0
SOURCE_GEOMETRY=1

OPERATION_THRESHOLD=0
OPERATION_SEGMENT=1
OPERATION_POSTPROCESS=2

SEGMENTATION_VOXELS=0
SEGMENTATION_AMINOACIDS=1
SEGMENTATION_DALTON=2
SEGMENTATION_AUTOMATIC=3

MORPHOLOGY_DILATION=0
MORPHOLOGY_EROSION=1
MORPHOLOGY_CLOSING=2
MORPHOLOGY_OPENING=3


class XmippProtCreateMask3D(ProtCreateMask3D, XmippGeometricalMask3D):
    """ Create a 3D mask.
    The mask can be created with a given geometrical shape (Sphere, Box, Cylinder...) or
    it can be obtained from operating on a 3d volume or a previuous mask.
    """
    _label = 'create 3d mask'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Mask generation')
        form.addParam('source', EnumParam, default=SOURCE_VOLUME,
                      choices=['Volume','Geometry'],
                      label='Mask source')
        # For volume sources
        isVolume = 'source==%d' % SOURCE_VOLUME
        form.addParam('inputVolume', PointerParam, pointerClass="Volume", condition=isVolume, 
                      label="Input volume",
                      help="Select the volume that will be used to create the mask")
        form.addParam('volumeOperation', EnumParam, default=OPERATION_THRESHOLD, condition=isVolume, 
                      choices=['Threshold','Segment','Only postprocess'],
                      label='Operation')
        #TODO: add wizard
        form.addParam('threshold', FloatParam, default=0.0,
                      condition='volumeOperation==%d and %s' % (OPERATION_THRESHOLD, isVolume),
                      label='Threshold')
        isSegmentation = 'volumeOperation==%d and %s' % (OPERATION_SEGMENT, isVolume)
        form.addParam('segmentationType', EnumParam, default=SEGMENTATION_DALTON, 
                      condition=isSegmentation,
                      label='Segmentation type',
                      choices=['Number of voxels','Number of aminoacids','Dalton mass','Automatic'])
        form.addParam('nvoxels', IntParam, 
                      condition='%s and segmentationType==%d' % (isSegmentation, SEGMENTATION_VOXELS),
                      label='Number of voxels')
        form.addParam('naminoacids', IntParam,
                      condition='%s and segmentationType==%d' % (isSegmentation, SEGMENTATION_AMINOACIDS),
                      label='Number of aminoacids')
        form.addParam('dalton', FloatParam, 
                      condition='%s and segmentationType==%d' % (isSegmentation, SEGMENTATION_DALTON),
                      label='Mass (Da)')
        
        # For geometrical sources
        form.addParam('samplingRate', FloatParam, default=1, 
                      condition='source==%d' % SOURCE_GEOMETRY, 
                      label="Sampling Rate (A/px)")
        XmippGeometricalMask3D.defineParams(self, form, 
                                            isGeometry='source==%d' % SOURCE_GEOMETRY, 
                                            addSize=True)

        # Postprocessing
        form.addSection(label='Postprocessing')
        form.addParam('doSmall', BooleanParam, default=False,
                      label='Remove small objects',
                      help="The input mask has to be binary")
        form.addParam('smallSize', IntParam, default=50,
                      label='Minimum size',condition="doSmall",
                      help="Connected components whose size is smaller than this number in voxels will be removed")
        form.addParam('doBig', BooleanParam, default=False,
                      label='Keep largest component',
                      help="The input mask has to be binary")
        form.addParam('doSymmetrize', BooleanParam, default=False,
                      label='Symmetrize mask')
        form.addParam('symmetry', StringParam, default='c1',
                      label='Symmetry group',condition="doSymmetrize",
                      help="See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry \n"
                           "for a description of the symmetry groups format. \n"
                           "If no symmetry is present, give c1")
        form.addParam('doMorphological', BooleanParam, default=False,
                      label='Apply morphological operation',
                      help="Dilation (dilate white region). \n"
                           "Erosion (erode white region). \n"
                           "Closing (Dilation+Erosion, removes black spots). \n"
                           "Opening (Erosion+Dilation, removes white spots). \n")
        form.addParam('morphologicalOperation', EnumParam, default=MORPHOLOGY_DILATION,
                      condition="doMorphological", 
                      choices=['dilation', 'erosion', 'closing', 'opening'],
                      label='Operation')
        form.addParam('elementSize', IntParam, default=1, condition="doMorphological",
                      expertLevel=LEVEL_ADVANCED,
                      label='Structural element size',                      
                      help="The larger this value, the more the effect will be noticed")
        form.addParam('doInvert', BooleanParam, default=False,
                      label='Invert the mask')
        form.addParam('doSmooth', BooleanParam, default=False,
                      label='Smooth borders',
                      help="Smoothing is performed by convolving the mask with a Gaussian.")
        form.addParam('sigmaConvolution', FloatParam, default=2, condition="doSmooth",
                      label='Gaussian sigma (px)',
                      help="The larger this value, the more the effect will be noticed")

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self.maskFile = self._getPath('mask.vol')
        
        if self.source == SOURCE_VOLUME:
            self._insertFunctionStep('createMaskFromVolumeStep')
        elif self.source == SOURCE_GEOMETRY:
            self.inputVolume.set(None)
            self._insertFunctionStep('createMaskFromGeometryStep')
        self._insertFunctionStep('postProcessMaskStep')
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions --------------------------------------------
    def createMaskFromVolumeStep(self):
        volume = self.inputVolume.get()
        fnVol = getImageLocation(volume)
        Ts = volume.getSamplingRate()
        
        if self.volumeOperation == OPERATION_THRESHOLD:
            self.runJob("xmipp_transform_threshold",
                        "-i %s -o %s --select below %f --substitute binarize" % (fnVol, self.maskFile,
                                                                                 self.threshold.get()))
        elif self.volumeOperation == OPERATION_SEGMENT:
            args="-i %s -o %s --method " % (fnVol, self.maskFile)
            if self.segmentationType == SEGMENTATION_VOXELS:
                args += "voxel_mass %d" % (self.nvoxels.get())
            elif self.segmentationType == SEGMENTATION_AMINOACIDS:
                args += "aa_mass %d %f" % (self.naminoacids.get(), Ts)
            elif self.segmentationType == SEGMENTATION_DALTON:
                args += "dalton_mass %d %f" % (self.dalton.get(), Ts)
            else:
                args += "otsu"
            self.runJob("xmipp_volume_segment", args)
        
        elif self.volumeOperation == OPERATION_POSTPROCESS:
            ImageHandler().convert(fnVol,self.maskFile)
            
        return [self.maskFile]
        
    def createMaskFromGeometryStep(self):
        # Create empty volume file with desired dimensions
        size = self.size.get()
        xmipp.createEmptyFile(self.maskFile, size, size, size)
        
        # Create the mask
        args = '-i %s ' % self.maskFile
        args += XmippGeometricalMask3D.argsForTransformMask(self,size)
        args += ' --create_mask %s' % self.maskFile
        self.runJob("xmipp_transform_mask", args)
        
        return [self.maskFile]
    
    def postProcessMaskStep(self):        
        if self.doSmall:
            self.runJob("xmipp_transform_morphology","-i %s --binaryOperation removeSmall %d"%(self.maskFile,self.smallSize.get()))
        
        if self.doBig:
            self.runJob("xmipp_transform_morphology","-i %s --binaryOperation keepBiggest"%self.maskFile)
        
        if self.doSymmetrize:
            if self.symmetry!='c1':
                self.runJob("xmipp_transform_symmetrize","-i %s --sym %s --dont_wrap"%(self.maskFile,self.symmetry.get()))
        
        if self.doMorphological:
            self.runJob("xmipp_transform_morphology","-i %s --binaryOperation %s --size %d"
                        %(self.maskFile,self.getEnumText('morphologicalOperation'),self.elementSize.get()))
        
        if self.doInvert:
            self.runJob("xmipp_image_operate","-i %s --mult -1"%self.maskFile)
            self.runJob("xmipp_image_operate","-i %s --plus  1"%self.maskFile)
        
        if self.doSmooth:
            self.runJob("xmipp_transform_filter","-i %s --fourier real_gaussian %f"%(self.maskFile,self.sigmaConvolution.get()))

    def createOutputStep(self):
        volMask = VolumeMask()
        volMask.setFileName(self.maskFile)
        
        if self.source==SOURCE_VOLUME:
            volMask.setSamplingRate(self.inputVolume.get().getSamplingRate())
        else:
            volMask.setSamplingRate(self.samplingRate.get())
        
        self._defineOutputs(outputMask=volMask)
        
        if self.source==SOURCE_VOLUME:
            self._defineSourceRelation(self.inputVolume, self.outputMask)
        
    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        messages = []      
        messages.append("*Mask creation*")
        if self.source==SOURCE_VOLUME:
            if self.volumeOperation==OPERATION_THRESHOLD:
                messages.append("   Thresholding %f"%self.threshold.get())
            elif self.volumeOperation==OPERATION_SEGMENT:
                if self.segmentationType==SEGMENTATION_AUTOMATIC:
                    messages.append("   Automatically segmented")
                else:
                    m="   Segmented to a mass of "
                    if self.segmentationType==SEGMENTATION_VOXELS:
                        m+="%d voxels"%(int(self.nvoxels.get()))
                    elif self.segmentationType==SEGMENTATION_AMINOACIDS:
                        m+="%d aminoacids"%(int(self.naminoacids.get()))
                    elif self.segmentationType==SEGMENTATION_DALTON:
                        m+="%d daltons"%(int(self.dalton.get()))
                    messages.append(m)
        elif self.source==SOURCE_GEOMETRY:
            size = self.size.get()
            messages.append("   Mask of size: %d x %d x %d"%(size,size,size))
            messages += XmippGeometricalMask3D.summary(self)

        messages.append("*Mask processing*")
        if self.doSmall:
            messages.append("   Removing components smaller than %d" % self.smallSize.get())
        if self.doBig:
            messages.append("   Keeping largest component")
        if self.doSymmetrize:
            messages.append("   Symmetrized %s" % self.symmetry.get())
        if self.doMorphological:
            messages.append("   Morphological operation: %s" % self.getEnumText('morphologicalOperation'))
        if self.doInvert:
            messages.append("   Inverted")
        if self.doSmooth:
            messages.append("   Smoothed (sigma=%f)"%self.sigmaConvolution.get())
        return messages

    def _citations(self):
        if (self.source == SOURCE_VOLUME and 
            self.volumeOperation == OPERATION_SEGMENT and 
            self.segmentationType==SEGMENTATION_AUTOMATIC):
            return ['Otsu1979']

    def _methods(self):
        messages = []      
        messages.append("*Mask creation*")
        
        if self.source == SOURCE_VOLUME:
            messages.append('We processed the volume %s.'%self.inputVolume.get().getNameId())

            if self.volumeOperation == OPERATION_THRESHOLD:
                messages.append("We thresholded it to a gray value of %f. "%self.threshold.get())
            elif self.volumeOperation == OPERATION_SEGMENT:
                if self.segmentationType == SEGMENTATION_AUTOMATIC:
                    messages.append("We automatically segmented it using Otsu's method [Otsu1979]")
                else:
                    m="We segmented it to a mass of "
                    if self.segmentationType == SEGMENTATION_VOXELS:
                        m+="%d voxels"%(int(self.nvoxels.get()))
                    elif self.segmentationType == SEGMENTATION_AMINOACIDS:
                        m+="%d aminoacids"%(int(self.naminoacids.get()))
                    elif self.segmentationType == SEGMENTATION_DALTON:
                        m+="%d daltons"%(int(self.dalton.get()))
                    messages.append(m)
        
        elif self.source == SOURCE_GEOMETRY:
            size=self.size.get()
            messages.append("We created a mask of size: %d x %d x %d voxels. "%(size,size,size))
            messages+=XmippGeometricalMask3D.methods(self)

        if self.doSmall:
            messages.append("We removed components smaller than %d voxels."%(self.smallSize.get()))
        if self.doBig:
            messages.append("We kept the largest component. ")
        if self.doSymmetrize:
            messages.append("We symmetrized it as %s. "%self.symmetry.get())
        if self.doMorphological:
            messages.append("Then, we applied a morphological operation, concisely, a %s. "%self.getEnumText('morphologicalOperation'))
        if self.doInvert:
            messages.append("We inverted the mask. ")
        if self.doSmooth:
            messages.append("And, we smoothed it (sigma=%f voxels)." % self.sigmaConvolution.get())
        if self.hasAttribute('outputMask'):
            messages.append('We refer to the output mask as %s.' % self.outputMask.getNameId())
        return messages
    