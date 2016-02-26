# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Josue Gomez Blanco     (jgomez@cnb.csic.es)
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
# import xmipp
from geometrical_mask import XmippGeometricalMask3D, XmippGeometricalMask2D
from protocol_process import XmippProcessParticles, XmippProcessVolumes
from pyworkflow.em.constants import *

from ..convert import getImageLocation

from ..constants import *


class XmippProtMask():
    """ This class implement the common features for applying a mask with Xmipp either SetOfParticles, Volume or SetOfVolumes objects.
    """
    
    def __init__(self, **args):
        self._program = "xmipp_transform_mask"
    
    #--------------------------- DEFINE param functions --------------------------------------------
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
        
        form.addParam('inputMask', PointerParam, pointerClass=self.MASK_CLASSNAME, 
                      condition='source==%d' % SOURCE_MASK,
                      label="Input mask")
        
        self.GEOMETRY_BASECLASS.defineParams(self, form, 
                                             isGeometry='source==%d' % SOURCE_GEOMETRY,
                                             addSize=False)
        
        form.addParam('fillType', EnumParam, 
                      choices=['value', 'min', 'max', 'avg'], 
                      condition='source==%d' % SOURCE_GEOMETRY,
                      default=MASK_FILL_VALUE,
                      label="Fill with ", display=EnumParam.DISPLAY_COMBO,
                      help='Select how are you going to fill the pixel values outside the mask. ')
        
        form.addParam('fillValue', IntParam, default=0, 
                      condition='fillType == %d and source==%d' % (MASK_FILL_VALUE,SOURCE_GEOMETRY),
                      label='Fill value',
                      help='Value to fill the pixel values outside the mask. ')
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertProcessStep(self):
        inputFn = self.inputFn

        if self.source == SOURCE_MASK:
            inputMaskFile = getImageLocation(self.inputMask.get())
            outputMaskFile = self._getTmpPath(self.MASK_FILE)
            self._insertFunctionStep('copyMaskFileStep', inputMaskFile, outputMaskFile)
        
        if self.fillType == MASK_FILL_VALUE:
            fillValue = self.fillValue.get()
        else:
            fillValue = self.getEnumText('fillType')
        
        if self.source == SOURCE_GEOMETRY:
            self._program = "xmipp_transform_mask"
            self._args += self._getGeometryCommand()
            self._args += " --substitute %(fillValue)s "
        elif self.source == SOURCE_MASK:
            self._program = "xmipp_image_operate"
            self._args += (" --mult %s" % outputMaskFile)
        else:
            raise Exception("Unrecognized mask type: %d" % self.source)
        
        self._insertFunctionStep("maskStep", self._args % locals())
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def copyMaskFileStep(self, inputMaskFile, outputMaskFile):
        """ Create a local copy of the mask.
        We use ImageHandler.convert instead of copyFile 
        because the mask could be inside an stack.
        """
        ImageHandler().convert(self.inputMask.get(), outputMaskFile)
    
    def maskStep(self, args):
        args += self._getExtraMaskArgs()
        self.runJob(self._program, args)            
        
    #--------------------------- UTILS functions ---------------------------------------------------
    def _getExtraMaskArgs(self):
        """ Return some extra arguments for the mask program.
        This function will varies when masking particles or volumes.
        """
        args = " -o %s --save_metadata_stack %s --keep_input_columns" % (self.outputStk, self.outputMd)
        return args

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self, geoClass):
        messages = []      
        messages.append("*Mask application*")
        if self.source.get() == SOURCE_GEOMETRY:
            messages.append(' Used geometrical mask:')
            messages += geoClass.summary(self)
        else:
            messages.append(' Used created mask: %s' % self.inputMask.get().getNameId())
            
        if self.fillType.get() == MASK_FILL_VALUE:
            messages.append(' Filled with %s value' % self.fillValue.get())
        else:
            messages.append(' Filled with %s value' % self.getEnumText('fillType'))
                                 
        return messages    

    def _methods(self, geoClass):
        messages = []      
        messages.append("*Mask application*")
        
        if self.source.get() == SOURCE_GEOMETRY:
            messages.append("We applied a geometrical mask:")
            messages+=geoClass.methods(self)
        else:
            messages.append('We applied a created mask: %s' % self.inputMask.get().getNameId())
            
        if self.fillType.get() == MASK_FILL_VALUE:
            messages.append('filled with %s value' % self.fillValue.get())
        else:
            messages.append('filled with %s value' % self.getEnumText('fillType'))

        return messages
    
    def _validateDimensions(self, inputSetName, inputMaskName, errorMsg):
        """ Validate that input set (either particles or volumens) have the same
        dimension than the input mask.
        """
        errors = []
        
        if self.source == SOURCE_MASK:
            px, py, _ = self.getAttributeValue(inputSetName).getDim()
            mx, my, _ = self.getAttributeValue(inputMaskName).getDim()
            if px != mx or py != my:
                errors.append(errorMsg)
                
        return errors        
        
        
class XmippProtMaskParticles(ProtMaskParticles, XmippProcessParticles, XmippProtMask, XmippGeometricalMask2D):
    """ Apply mask to a set of particles """
    _label = 'apply 2d mask'
    
    MASK_FILE = 'mask.spi'
    MASK_CLASSNAME = 'Mask'
    GEOMETRY_BASECLASS = XmippGeometricalMask2D
    
    def __init__(self, **kwargs):
        ProtMaskParticles.__init__(self, **kwargs)
        XmippProcessParticles.__init__(self, **kwargs)
        XmippProtMask.__init__(self, **kwargs)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        XmippProtMask._defineProcessParams(self, form)
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _getGeometryCommand(self):
        Xdim = self.inputParticles.get().getDimensions()[0]
        self.ndim = self.inputParticles.get().getSize()
        args = XmippGeometricalMask2D.argsForTransformMask(self, Xdim)
        return args

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        messages = []      
        messages += XmippProtMask._summary(self, XmippGeometricalMask2D)

        return messages

    def _methods(self):
        messages = []      
        messages += XmippProtMask._methods(self, XmippGeometricalMask2D)

        return messages
    
    def _validate(self):
        return XmippProtMask._validateDimensions(self, 
                                                 'inputParticles', 'inputMask',
                                                 'Input particles and mask should '
                                                 'have same dimensions.')
    
    
class XmippProtMaskVolumes(ProtMaskVolumes, XmippProcessVolumes, XmippProtMask, XmippGeometricalMask3D):
    """ Apply mask to a volume """
    _label = 'apply 3d mask'
    
    MASK_FILE = 'mask.vol'
    MASK_CLASSNAME = 'VolumeMask'
    GEOMETRY_BASECLASS = XmippGeometricalMask3D
    
    def __init__(self, **kwargs):
        ProtMaskVolumes.__init__(self, **kwargs)
        XmippProcessVolumes.__init__(self, **kwargs)
        XmippProtMask.__init__(self, **kwargs)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        XmippProtMask._defineProcessParams(self, form)
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _getGeometryCommand(self):
        if self._isSingleInput():            
            Xdim = self.inputVolumes.get().getDim()[0]
        else:
            Xdim = self.inputVolumes.get().getDimensions()[0]
        args = XmippGeometricalMask3D.argsForTransformMask(self, Xdim)
        return args

    def _getExtraMaskArgs(self):
        if self._isSingleInput():
            return " -o %s" % self.outputStk
        else:
            return XmippProtMask._getExtraMaskArgs(self)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        messages = []      
        messages += XmippProtMask._summary(self, XmippGeometricalMask3D)

        return messages
    
    def _methods(self):
        messages = []      
        messages += XmippProtMask._methods(self, XmippGeometricalMask3D)

        return messages
    
    def _validate(self):
        return XmippProtMask._validateDimensions(self, 
                                                 'inputVolumes', 'inputMask',
                                                 'Input volumes and mask should '
                                                 'have same dimensions.')
        