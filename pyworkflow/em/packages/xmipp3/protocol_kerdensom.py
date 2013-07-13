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
        self._addParams()
        
    def _addParams(self):
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
                      'See [http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM]')
        self.addParam('SomReg1', IntParam, default=200, expertLevel=LEVEL_ADVANCED,
                      label='Final regularization factor:')
        self.addParam('SomSteps', IntParam, default=5, expertLevel=LEVEL_ADVANCED,
                      label='Regularization steps:',
                      help='Number of steps to lower the regularization factor')
        self.addParam('extraParams', StringParam,
                      label="Additional kerdenSOM parameters:", 
                      help='For a complete description'
                      'See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM')
        
        
class XmippProtKerdensom(ProtClassify, XmippProtocol):
    """ Protocol to align a set of particles. """
    _definition = XmippDefKerdensom()
    _label = 'Xmipp KerDenSom'

    def _defineSteps(self):
        self._prepareDefinition()       
        self._insertSteps()
    
    
    def _prepareDefinition(self):
         # Convert input images if necessary
        self.inputImgs = self.inputImages.get()        
        imgsFn = self._insertConvertStep('inputImgs', XmippSetOfParticles,
                                         self._getPath('input_images.xmd'))
        mask = self.Mask.get()
        self._params = {'oroot': self._getExtraPath("kerdensom"),
                        'imgsFn': imgsFn,
                        'mask': mask,
                        'SomXdim': self.SomXdim.get(),
                        'SomYdim': self.SomYdim.get(),
                        'SomReg0': self.SomReg0.get(),
                        'SomReg1': self.SomReg1.get(),
                        'SomSteps': self.SomSteps.get(),
                        'extraParams': self.extraParams.get(),
                        'vectors': self._getExtraPath("vectors.xmd"),
                        'classes': self._getExtraPath("classes.stk"),
                        'kvectors': self._getExtraPath("kerdensom_vectors.xmd"),
                        'kclasses': self._getExtraPath("kerdensom_classes.xmd")
                       }    
        
    def _insertSteps(self):
        self._insertImgToVector()
        self._insertKerdensom()
        self._insertVectorToImg()
        self._insertFunctionStep('rewriteClassBlock')
        self._insertFunctionStep('createOutput')
        
    def _insertImgToVector(self):
        """ Insert runJob for convert into a vector Md """
        args = ' -i %(imgsFn)s -o %(vectors)s '
        if self.useMask:
            args += ' --mask binary_file %(mask)s'
        
        self._insertRunJobStep("xmipp_image_vectorize", args % self._params)

    def _insertKerdensom(self):
        args = '-i %(vectors)s --oroot %(oroot)s --xdim %(SomXdim)d --ydim %(SomYdim)d' + \
               ' --deterministic_annealing %(SomSteps)f %(SomReg0)f %(SomReg1)f %(extraParams)s'
        self._insertRunJobStep("xmipp_classify_kerdensom", args % self._params)
#        deleteFiles([self._getExtraPath("vectors.xmd"),self._getExtraPath("vectors.vec")], True)
   
    def _insertVectorToImg(self):
        args = ' -i %(kvectors)s -o %(classes)s' 
        if self.useMask:
            args += ' --mask binary_file %(mask)s'
        self._insertRunJobStep("xmipp_image_vectorize", args % self._params)
#        deleteFiles([self._getExtraPath("kerdensom_vectors.xmd"),self._getExtraPath("kerdensom_vectors.vec")], True)

    def rewriteClassBlock(self):
        fnClass = "classes@%(kclasses)s" % self._params
        fnClassStack = self._params['classes']
        md = xmipp.MetaData(fnClass)
        # TODO: Check if following is necessary
        counter = 1
        for id in md:
            md.setValue(xmipp.MDL_IMAGE,"%06d@%s"%(counter,fnClassStack),id)
            counter += 1
        md.write(fnClass, xmipp.MD_APPEND)
        #createLink(Log, fnClassMetadata, self.getPath("classes.xmd"))
        #createLink(Log, self._getExtraPath("kerdensom_images.xmd"), self.getPath("images.xmd"))
    
    def createOutput(self):
        """ Store the kenserdom object 
        as result of the protocol. 
        """
        classification = XmippClassification2D(self._params['kclasses'])
        self._defineOutputs(outputClassification=classification)

    def _validate(self):
        errors = []
        mask = self.Mask.get()
        if self.SomReg0 < self.SomReg1:
            errors.append("Regularization must decrease over iterations:")
            errors.append("    Initial regularization must be larger than final")
        if self.useMask:
            if len(mask) > 0:
                if not exists(mask):
                    errors.append("Cannot find the file " + mask)
            else:
                errors.append("Please, enter a mask file")
        return errors
        
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputClassification'):
            summary.append("Output classification not ready yet.")
        else:
            summary.append("Input Images: %s" % self.inputImages.get().getNameId())
            summary.append("Output Classified Images: %s" % self.outputClassification.get())
        return summary
