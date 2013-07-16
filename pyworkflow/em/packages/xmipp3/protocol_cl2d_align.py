# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-package contains wrapper around align2d Xmipp program
"""

from os.path import join, dirname, exists
from pyworkflow.em import *  
import xmipp
from data import *
from xmipp3 import XmippProtocol
from glob import glob

class XmippDefCL2DAlign(Form):
    """Create the definition of parameters for
    the XmippProtAlign2d protocol"""
    def __init__(self):
        Form.__init__(self)

        self.addSection(label='Input')
        self.addParam('inputImages', PointerParam, label="Input images", important=True, 
                      pointerClass='SetOfParticles',
                      help='Select the input images from the project.'
                           'It should be a SetOfImages class')
        
        self.addParam('useReferenceImage', BooleanParam, default=True,
                      label='Use a Reference Image ?', 
                      help='If you set to <Yes>, you should provide a reference image'
                           'If <No>, the default generation is done by averaging'
                           'subsets of the input images.')
        self.addParam('ReferenceImage', StringParam , condition='useReferenceImage',
                      label="Reference image(s)", 
                      help='Image that will serve as class reference')
        
        self.addParam('maximumShift', IntParam, default=10,
                      label='Maximum shift:',
                      help='in pixels')
        self.addParam('numberOfIterations', IntParam, default=10, expertLevel=LEVEL_ADVANCED,
                      label='Number of iterations:',
                      help='Maximum number of iterations')
        self.addParallelSection(threads=1, mpi=3)
        
        
class XmippProtCL2DAlign(ProtAlign, XmippProtocol):
    """ Protocol to align a set of particles. """
    _definition = XmippDefCL2DAlign()
    _label = 'Xmipp CL2D Align'

    def _defineSteps(self):
        """ Mainly prepare the command line for call cl2d align program"""
        
        # Convert input images if necessary
        self.inputImgs = self.inputImages.get()        
        imgsFn = self._insertConvertStep('inputImgs', XmippSetOfImages,
                                         self._getPath('input_images.xmd'))
        # Prepare arguments to call program: xmipp_classify_CL2D
        self._params = {'imgsFn': imgsFn, 
                        'extraDir': self._getExtraPath(),
                        'maxshift': self.maximumShift.get(),
                        'iter': self.numberOfIterations.get(),
                        }
        args = '-i %(imgsFn)s --odir %(extraDir)s --nref 1 --iter %(iter)d --maxShift %(maxshift)d'
        if self.useReferenceImage:
            args += " --ref0 " + self.ReferenceImage.get()
        else:
            args += " --nref0 1"
            
        self._defineClassifySteps("xmipp_classify_CL2D", args)
              
    def _defineClassifySteps(self, program, args, subset=''):
        self._insertRunJobStep(program, args % self._params)
        self._insertFunctionStep('createOutput')        
        
    def createOutput(self):
        """ Store the XmippClassification2D object 
        as result of the protocol. 
        """
        lastMdFn = self._getExtraPath("images.xmd")
        
        imgs = XmippSetOfParticles(lastMdFn)
        imgs.setAligment(True)
        self._defineOutputs(outputParticles=imgs)

    def validate(self):
        errors = []
        refImage = self.ReferenceImage.get()
        if self.useReferenceImage:
            if len(refImage) > 0:
                if not exists(refImage):
                    errors.append("Cannot find the file " + refImage)
            else:
                errors.append("Please, enter an Image file")
        return errors
        
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output alignment not ready yet.")
        else:
            summary.append("Input Images: %s" % self.inputImages.get().getNameId())
            summary.append("Output Aligned Images: %s" % self.outputParticles.get())
        return summary
