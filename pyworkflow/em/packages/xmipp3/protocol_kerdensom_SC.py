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
This sub-package contains wrapper around kendersom Xmipp program
"""

from os.path import join, dirname, exists
from pyworkflow.em import *  
import xmipp
from data import *
from xmipp3 import XmippProtocol
from glob import glob

class XmippDefKerdensom(Form):
    """Create the definition of parameters for
    the XmippProtAlign2d protocol"""
    def __init__(self):
        Form.__init__(self)

        self.addSection(label='Input')
        self.addParam('inputImages', PointerParam, label="Input images", important=True, 
                      pointerClass='SetOfParticles',
                      help='Select the input images from the project.'
                           'It should be a SetOfParticles class')
        self.addParam('useMask', BooleanParam, default=False,
                      label='Use a Mask ?', 
                      help='If you set to <Yes>, you should provide a mask')
        self.addParam('Mask', StringParam , condition='useMask',
                      label="Mask", 
                      help='Mask image will serve to enhance the classification')
        self.addParam('SomXdim', IntParam, default=7,
                      label='X-dimension of the map:')
        self.addParam('SomYdim', IntParam, default=7,
                      label='Y-dimension of the map:')
        self.addParam('SomReg0', IntParam, default=1000, expertLevel=LEVEL_ADVANCED,
                      label='Initial regularization factor', 
                      help='The kerdenSOM algorithm anneals from an initial high regularization factor'
                      'to a final lower one, in a user-defined number of steps.'
                      'If the output map is too smooth, lower the regularization factors'
                      'If the output map is not organized, higher the regularization factors'
                      'See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM')
        self.addParam('SomReg1', IntParam, default=200, expertLevel=LEVEL_ADVANCED,
                      label='Final regularization factor:')
        self.addParam('SomSteps', IntParam, default=5, expertLevel=LEVEL_ADVANCED,
                      label='Regularization steps:',
                      help='Number of steps to lower the regularization factor')
        self.addParam('KerdensomExtraCommand', StringParam,
                      label="Additional kerdenSOM parameters:", 
                      help='For a complete description'
                      'See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM')
        
        
class XmippProtKerdensom(ProtClassify, XmippProtocol):
    """ Protocol to align a set of particles. """
    _definition = XmippDefKerdensom()
    _label = 'Xmipp KerDenSom'

    def _defineSteps(self):
        # Convert input images if necessary
        self.inputImgs = self.inputImages.get()        
        imgsFn = self._insertConvertStep('inputImgs', XmippSetOfParticles,
                                         self._getPath('input_images.xmd'))
        self._params = {'extraDir': self._getExtraPath(),
                        'mask': self.Mask.get(), 
                        'SomXdim': self.SomXdim.get(),
                        'SomYdim': self.SomYdim.get(),
                        'SomReg0': self.SomReg0.get(),
                        'SomReg1': self.SomReg1.get(),
                        'SomSteps': self.SomSteps.get(),
                        'KerdensomExtraCommand': self.KerdensomExtraCommand.get()
                       }
        args1 = '-i ' + self._getExtraPath("vectors.xmd") + ' --oroot ' + self._getExtraPath("kerdensom") + ' --xdim %(SomXdim)d --ydim %(SomYdim)d' + \
                ' --deterministic_annealing %(SomSteps)f %(SomReg0)f %(SomReg1)f '
        if self.KerdensomExtraCommand != '':
            args1 += '%(KerdensomExtraCommand)s'
       
        self._insertFunctionStep('img2vector', imgsFn)
        self._insertFunctionStep('kerdensom', args1 % self._params)
        self._insertFunctionStep('vector2img')
        self._insertFunctionStep('rewriteClassBlock')
        self._insertFunctionStep('createOutput')
        
    def img2vector(self, imgsFn):
        """convert into a vector Md """
        mask = self._params['mask']
        args = ' -i ' + imgsFn + ' -o ' + self._getExtraPath("vectors.xmd")
        if mask != '':
            args += ' --mask binary_file ' + mask
        self.runJob(None,"xmipp_image_vectorize", args)

    def kerdensom(self, args1):
        self.runJob(None,"xmipp_classify_kerdensom", args1)
#        deleteFiles([self._getExtraPath("vectors.xmd"),self._getExtraPath("vectors.vec")], True)
   
    def vector2img(self):
        mask = self._params['mask']
        args = ' -i ' + self._getExtraPath("kerdensom_vectors.xmd") + ' -o ' + self._getExtraPath("classes.stk")
        if mask != '':
            args += ' --mask binary_file ' + mask
        self.runJob(None,"xmipp_image_vectorize", args)
#        deleteFiles([self._getExtraPath("kerdensom_vectors.xmd"),self._getExtraPath("kerdensom_vectors.vec")], True)

def rewriteClassBlock(self):
    fnClassMetadata = self._getExtraPath("kerdensom_classes.xmd")
    fnClass = "classes@%s" %fnClassMetadata
    fnClassStack = self._getExtraPath("classes.stk")
    mD = MetaData(fnClass)
    counter = 1
    for id in mD:
        mD.setValue(MDL_IMAGE,"%06d@%s"%(counter,fnClassStack),id)
        counter += 1
    mD.write(fnClass, xmipp.MD_APPEND)
    createLink(Log, fnClassMetadata, self.getPath("classes.xmd"))
    createLink(Log, self._getExtraPath("kerdensom_images.xmd"), self.getPath("images.xmd"))
    
    def createOutput(self):
        """ Store the kenserdom object 
        as result of the protocol. 
        """
        classification = XmippClassification2D(self.oroot + 'classes.xmd')
        self._defineOutputs(outputClassification=classification)

    def validate(self):
        errors = []
        if self.SomReg0 < self.SomReg1:
            errors.append("Regularization must decrease over iterations:")
            errors.append("    Initial regularization must be larger than final")
        return errors
        
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputClassification'):
            summary.append("Output alignment not ready yet.")
        else:
            summary.append("Input Images: %s" % self.inputImages.get().getNameId())
            summary.append("Output Aligned Images: %s" % self.outputClassification.get())
        return summary
