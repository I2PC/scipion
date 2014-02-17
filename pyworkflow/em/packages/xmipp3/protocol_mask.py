# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
This sub-package contains protocol for masks operations
"""

from pyworkflow.em import *  
from pyworkflow.utils import *  
import xmipp
from geometrical_mask import XmippGeometricalMask3D, XmippGeometricalMask2D
from protocol_process import XmippProcessParticles, XmippProcessVolumes
from convert import createXmippInputImages, readSetOfParticles

from pyworkflow.em.constants import *
from constants import *


class XmippProtMask():
    """ This class implement a protocol for applying a mask with Xmipp.
    """
    
    def __init__(self, **args):
        self._program = "xmipp_transform_mask"

    def _defineProcessParams(self, form):
        """ Add common mask parameters that can be used
        in several protocols definitions.
        Params:
            form: the Definition instance.
        """
        
        form.addParam('source', EnumParam,
                      label='Mask source',
                      default=SOURCE_GEOMETRY, choices=['Geometry','Created mask'], 
                      help='Select which type of mask do you want to apply. \n ')
        
        self._defineProtParams(form)
        
        form.addParam('fillType', EnumParam, 
                      choices=['value', 'min', 'max', 'avg'], 
                      default=MASK_FILL_VALUE,
                      label="Fill with ", display=EnumParam.DISPLAY_COMBO,
                      help='Select how are you going to fill the pixel values outside the mask. ')
        
        form.addParam('fillValue', IntParam, default=0, 
                      condition='fillType == %d' % MASK_FILL_VALUE,
                      label='Fill value',
                      help='Value to fill the pixel values outside the mask. ')
    
    def _getCommand(self, inputFn):
        fillStr = self.getEnumText('fillType')
#         maskRadius = self.maskRadius.get()
#         maskBand = maskRadius + self.maskOuterRadius.get()
#         
#         # Weird logic of Xmipp, minus to define the mask
#         maskRadius *= -1
#         maskBand *= -1
        
        if self.fillType == MASK_FILL_VALUE:
            fillStr = str(self.fillValue.get())
        
        self._args += " --substitute %(fillStr)s "
        
        if self.source == SOURCE_GEOMETRY:
            self._args += self._getGeometryCommand()
        elif self.source == SOURCE_MASK:
            self._args += "--mask binary_file %s" % self.maskFn
        else:
            raise Exception("Unrecognized mask type: %d" % self.source.get())

        return self._args % locals()


class XmippProtMaskParticles(ProtMaskParticles, XmippProtMask, XmippProcessParticles, XmippGeometricalMask2D):
    """ Apply some filter to SetOfParticles """
    _label = 'mask particles'
    
    def __init__(self, **args):
        ProtMaskParticles.__init__(self)
        XmippProtMask.__init__(self, **args)
        XmippProcessParticles.__init__(self, **args)
#         XmippGeometricalMask2D.__init__(self, **args)
        
    def _defineProtParams(self, form):
        form.addParam('inputMask', PointerParam, pointerClass="Mask", label="Input mask",condition='source==%d'%SOURCE_MASK)
        XmippGeometricalMask2D.defineParams(self, form, isGeometry='source==%d' % SOURCE_GEOMETRY, addSize=False)
    
    def _defineFilenames(self):
        XmippProcessParticles._defineFilenames(self)
        self.maskFn = self._getPath('mask.spi')
    
    def copyFileStep(self, *args):
        """ Convert the input mask to file. """
        ImageHandler().convert(self.inputMask.get().getLocation(), (None, self.maskFn))
    
    def _insertProcessStep(self, inputFn, outputFn, outputMd):
        if self.source == SOURCE_MASK:
            self._insertFunctionStep('copyFileStep', self.inputMask.get().getLocation())
        XmippProcessParticles._insertProcessStep(self, inputFn, outputFn, outputMd)
    
    def _getGeometryCommand(self):
        (Xdim, _, _, self.ndim) = self.inputParticles.get().getDimensions()
        args = XmippGeometricalMask2D.argsForTransformMask(self,Xdim)
        return args


class XmippProtMaskVolumes(ProtMaskVolumes, XmippProtMask, XmippProcessVolumes, XmippGeometricalMask3D):
    """ Apply mask to volume or SetOfVolumes """
    _label = 'apply mask'
    
    def __init__(self, **args):
        ProtMaskVolumes.__init__(self)
        XmippProtMask.__init__(self, **args)
        XmippProcessVolumes.__init__(self, **args)
    
    def _defineProtParams(self, form):
        form.addParam('inputMask', PointerParam, pointerClass="VolumeMask", label="Input mask",condition='source==%d'%SOURCE_MASK)
        XmippGeometricalMask3D.defineParams(self, form, isGeometry='source==%d'%SOURCE_GEOMETRY, addSize=False)
    
    def _defineFilenames(self):
        XmippProcessVolumes._defineFilenames(self)
        self.maskFn = self._getPath('mask.spi')
    
    def _insertProcessStep(self, inputFn, outputFn, outputMd):
        if self.source == SOURCE_MASK:
            self._insertFunctionStep('copyFileStep', self.inputMask.get().getLocation())
        XmippProcessVolumes._insertProcessStep(self, inputFn, outputFn, outputMd)
    
    def _getGeometryCommand(self):
        if isinstance(self.inputVolumes.get(), Volume):
            Xdim, _, _, _ = self.inputVolumes.get().getDim()
        else:
            Xdim, _, _, _ = self.inputVolumes.get().getDimensions()
        args = XmippGeometricalMask3D.argsForTransformMask(self,Xdim)
        return args

