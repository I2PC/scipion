# **************************************************************************
# *
# * Authors:     Laura del Cano (ldelcano@cnb.csic.es)
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

SOURCE_PARTICLE=0
SOURCE_GEOMETRY=1
SOURCE_MASK=2


class XmippProtCreateMask2D(ProtCreateMask2D, XmippGeometricalMask2D):
    """ Create a 2D mask.
    The mask can be created with a given geometrical shape (Circle, Rectangle, Crown...) or
    it can be obtained from operating on a 2d image or a previuous mask.
    """
    _label = 'create 2d mask'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Mask generation')
        
        # For geometrical sources
        form.addParam('samplingRate', FloatParam, default=1, 
                      label="Sampling Rate (A/px)")
        
        XmippGeometricalMask2D.defineParams(self, form, 
                                            isGeometry=True, 
                                            addSize=True)


    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self.maskFile = self._getPath('mask.xmp')
    
        self._insertFunctionStep('createMaskFromGeometryStep')

        self._insertFunctionStep('postProcessMaskStep')
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions --------------------------------------------
        
    def createMaskFromGeometryStep(self):
        # Create empty volume file with desired dimensions
        size = self.size.get()
        xmipp.createEmptyFile(self.maskFile, size, size)
        
        # Create the mask
        args = '-i %s ' % self.maskFile
        args += XmippGeometricalMask2D.argsForTransformMask(self,size)
        args += ' --create_mask %s' % self.maskFile
        self.runJob("xmipp_transform_mask", args)
        
        return [self.maskFile]


    def postProcessMaskStep(self):
        pass

    def createOutputStep(self):
        mask = Mask()
        mask.setFileName(self.maskFile)
        
        mask.setSamplingRate(self.samplingRate.get())
        
        self._defineOutputs(outputMask=mask)
        
    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        messages = []      
        messages.append("*Mask creation*")
        size = self.size.get()
        messages.append("   Mask of size: %d x %d"%(size,size))
        messages += XmippGeometricalMask2D.summary(self)

        return messages

    def _citations(self):
        pass

    def _methods(self):
        messages = []      
        messages.append("*Mask creation*")
        
        size=self.size.get()
        messages.append("We created a mask of size: %d x %d pixels. "%(size,size))
        messages+=XmippGeometricalMask2D.methods(self)

        if self.hasAttribute('outputMask'):
            messages.append('Output mask: %s.' % self.getObjectTag('outputMask'))
        return messages
    