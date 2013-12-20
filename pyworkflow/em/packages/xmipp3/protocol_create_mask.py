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
# *  e-mail address 'xmipp@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-package contains protocols for creating 3D masks.
"""

from pyworkflow.em import *  
from constants import *
from geometrical_mask import *

SOURCE_VOLUME=0
SOURCE_GEOMETRY=1
SOURCE_MASK=2

OPERATION_THRESHOLD=0
OPERATION_SEGMENT=1

SEGMENTATION_VOXELS=0
SEGMENTATION_AMINOACIDS=1
SEGMENTATION_DALTON=2
SEGMENTATION_AUTOMATIC=3

MORPHOLOGY_DILATION=0
MORPHOLOGY_EROSION=1
MORPHOLOGY_CLOSING=2
MORPHOLOGY_OPENING=3

class XmippProtCreateMask3D(ProtCreateMask3D,XmippGeometricalMask):
    """ Create a 3D mask from a geometrical description (Sphere, Box, Cylinder...), from a volume or from another class """
    _label = 'create mask'
    
    def _defineParams(self, form):
        form.addSection(label='Mask generation')
        form.addParam('source', EnumParam, label='Mask source', default=SOURCE_VOLUME, choices=['Volume','Geometry','Another mask'])
        
        # For volume sources
        isVolume = 'source==%d'%SOURCE_VOLUME
        form.addParam('volume', PointerParam, pointerClass="SetOfVolumes", label="Input volume",
                      condition=isVolume, help="Volume that will serve as basis for the mask")
        form.addParam('volumeOperation',EnumParam,label='Operation',condition=isVolume,
                      default=OPERATION_THRESHOLD,choices=['Threshold','Segment'])
        form.addParam('threshold',FloatParam,label='Threshold',condition='%s and volumeOperation==%d'%(isVolume,OPERATION_THRESHOLD))
        isSegmentation='%s and volumeOperation==%d'%(isVolume,OPERATION_SEGMENT)
        form.addParam('segmentationType',EnumParam,label='Segmentation type',condition=isSegmentation,
                      default=SEGMENTATION_DALTON,choices=['Number of voxels','Number of aminoacids','Dalton mass','Automatic'])
        form.addParam('nvoxels',IntParam,label='Number of voxels',
                      condition='%s and segmentationType==%d'%(isSegmentation,SEGMENTATION_VOXELS))
        form.addParam('naminoacids',IntParam,label='Number of aminoacids',
                      condition='%s and segmentationType==%d'%(isSegmentation,SEGMENTATION_AMINOACIDS))
        form.addParam('dalton',FloatParam,label='Mass (Da)',
                      condition='%s and segmentationType==%d'%(isSegmentation,SEGMENTATION_DALTON))
        
        # For geometrical sources
        XmippGeometricalMask.defineParams(self, form, 'source==%d'%SOURCE_GEOMETRY)

        # For another mask
        form.addParam('inputMask', PointerParam, pointerClass="VolumeMask", label="Input mask",condition='source==%d'%SOURCE_MASK)

        # Postprocessing
        form.addSection(label='Postprocessing')
        form.addParam('doSmall',BooleanParam,default=False,label='Remove small objects',help="The input mask has to be binary")
        form.addParam('smallSize',IntParam,default=50,label='Minimum size',condition="doSmall",
                      help="Connected components whose size is smaller than this number in voxels will be removed")
        form.addParam('doBig',BooleanParam,default=False,label='Keep largest component',help="The input mask has to be binary")
        form.addParam('doSymmetrize',BooleanParam,default=False,label='Symmetrize mask')
        form.addParam('symmetry',StringParam,default='c1',label='Symmetry group',condition="doSymmetrize",
                      help="See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry for a description of the symmetry groups format."
                           "If no symmetry is present, give c1")
        form.addParam('doMorphological',BooleanParam,default=False,label='Apply morphological operation',
                      help="Dilation (dilate white region)."
                           "Erosion (erode white region)."
                           "Closing (Dilation+Erosion, removes black spots)."
                           "Opening (Erosion+Dilation, removes white spots).")
        form.addParam('morphologicalOperation', EnumParam, label='Operation', default=MORPHOLOGY_DILATION,condition="doMorphological",
                      choices=['dilation', 'erosion', 'closing', 'opening'])
        form.addParam('elementSize',IntParam,default=1,label='Structural element size',condition="doMorphological",
                      expertLevel=LEVEL_ADVANCED,help="The larger this value, the more the effect will be noticed")
        form.addParam('doInvert',BooleanParam,default=False,label='Invert the mask')
        form.addParam('doSmooth',BooleanParam,default=False,label='Smooth borders',
                      help="Smoothing is performed by convolving the mask with a Gaussian.")
        form.addParam('sigmaConvolution',FloatParam,default=2,label='Gaussian sigma (px)',condition="doSmooth",
                      help="The larger this value, the more the effect will be noticed")

    def _defineSteps(self):
        self.maskFile = self._getPath('mask.vol')
        if self.source==SOURCE_VOLUME:
            self.inputMask.set(None)
            self._insertFunctionStep('createMaskFromVolume')
        elif self.source==SOURCE_GEOMETRY:
            self.inputMask.set(None)
            self.volume.set(None)
            self._insertFunctionStep('createMaskFromGeometry')
        else:
            self.volume.set(None)
            self._insertFunctionStep('createMaskFromAnotherMask')
        self._insertFunctionStep('postProcessMask')
        self._insertFunctionStep('createOutput')
    
    def createMaskFromVolume(self):
        volume=self.volume.get().getFirstItem()
        if self.volumeOperation==OPERATION_THRESHOLD:
            self.runJob(None,"xmipp_transform_threshold",
                        "-i %s -o %s --select below %f --substitute binarize"%(volume.getFileName(),self.maskFile,
                                                                               self.threshold.get()))
        elif self.volumeOperation==OPERATION_SEGMENT:
            args="-i %s -o %s --method "%(volume.getFileName(),self.maskFile)
            Ts=volume.getSamplingRate()
            if self.segmentationType==SEGMENTATION_VOXELS:
                args+="voxel_mass %d"%(self.nvoxels.get())
            elif self.segmentationType==SEGMENTATION_AMINOACIDS:
                args+="aa_mass %d %f"%(self.naminoacids.get(),Ts)
            elif self.segmentationType==SEGMENTATION_DALTON:
                args+="dalton_mass %d %f"%(self.dalton.get(),Ts)
            else:
                args+="otsu"
            self.runJob(None,"xmipp_volume_segment",args)
        return [self.maskFile]
        
    def createMaskFromGeometry(self):
        # Create empty volume file with desired dimensions
        size=self.size.get()
        xmipp.createEmptyFile(self.maskFile, size, size, size)
        
        # Create the mask
        args = '-i %s ' % self.maskFile
        args+=XmippGeometricalMask.argsForTransformMask(self)
        args += ' --create_mask %s' % self.maskFile
        self.runJob(None, "xmipp_transform_mask", args)
        return [self.maskFile]
    
    def createMaskFromAnotherMask(self):
        ImageHandler().convert(self.inputMask.get().getLocation(), (0, self.maskFile))

    def postProcessMask(self):
        if self.doSmall.get():
            self.runJob(None,"xmipp_transform_morphology","-i %s --binaryOperation removeSmall %d"%(self.maskFile,self.smallSize.get()))
        if self.doBig.get():
            self.runJob(None,"xmipp_transform_morphology","-i %s --binaryOperation keepBiggest"%self.maskFile)
        if self.doSymmetrize.get():
            if self.symmetry.get()!='c1':
                self.runJob(None,"xmipp_transform_symmetrize","-i %s --sym %s --dont_wrap"%(self.maskFile,self.symmetry.get()))
        if self.doMorphological.get():
            self.runJob(None,"xmipp_transform_morphology","-i %s --binaryOperation %s --size %d"
                        %(self.maskFile,self.getEnumText('morphologicalOperation'),self.elementSize.get()))
        if self.doInvert.get():
            self.runJob(None,"xmipp_image_operate","-i %s --mult -1"%self.maskFile)
            self.runJob(None,"xmipp_image_operate","-i %s --plus  1"%self.maskFile)
        if self.doSmooth.get():
            self.runJob(None,"xmipp_transform_filter","-i %s --fourier real_gaussian %f"%(self.maskFile,self.sigmaConvolution.get()))

    def createOutput(self):
        volMask = VolumeMask()
        volMask.setFileName(self.maskFile)
        if self.source==SOURCE_VOLUME:
            volMask.setSamplingRate(self.volume.get().getSamplingRate())
        if self.source==SOURCE_MASK:
            volMask.setSamplingRate(self.inputMask.get().getSamplingRate())
        self._defineOutputs(outputMask=volMask)
        
    def _summary(self):
        messages = []      
        messages.append("*Mask creation*")
        if self.source==SOURCE_MASK:
            messages.append("   From another mask")
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
            size=self.size.get()
            messages.append("Mask of size: %d x %d x %d"%(size,size,size))
            messages+=XmippGeometricalMask.summary(self)

        messages.append("*Mask processing*")
        if self.doSmall.get():
            messages.append("   Removing components smaller than %d"%(self.smallSize.get()))
        if self.doBig.get():
            messages.append("   Keeping largest component")
        if self.doSymmetrize.get():
            messages.append("   Symmetrized %s"%self.symmetry.get())
        if self.doMorphological.get():
            messages.append("   Morphological operation: %s"%self.getEnumText('morphologicalOperation'))
        if self.doInvert.get():
            messages.append("   Inverted")
        if self.doSmooth.get():
            messages.append("   Smoothed (sigma=%f)"%self.sigmaConvolution.get())
        return messages
    