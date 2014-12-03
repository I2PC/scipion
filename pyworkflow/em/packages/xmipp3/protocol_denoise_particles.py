# **************************************************************************
# *
# * Authors:     I. Foche (ifoche@cnb.csic.es)
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
This sub-package contains wrapper around Screen Particles Xmipp program
"""

from pyworkflow.object import String
from pyworkflow.protocol.params import (EnumParam, IntParam, Positive, Range,
                                        LEVEL_EXPERT, FloatParam, PointerParam)
from pyworkflow.em.protocol import ProtProcessParticles
from pyworkflow.utils.path import copyFile, replaceBaseExt

from convert import writeSetOfParticles, readSetOfParticles



        
        
class XmippProtDenoiseParticles(ProtProcessParticles):
    """ Remove particles noise using different Xmipp programs """
    _label = 'denoise particles'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        # First we customize the inputParticles param to fit our needs in this protocol
        form.getParam('inputParticles').pointerCondition = 'hasAlignment'
        form.getParam('inputParticles').help = String('Points to your input images. It is'' important that the images have alignment information with '
                                                                    ' important that the images have alignment information with '
                                                                    'respect to the chosen set of classes. This is the standard situation '
                                                                    'after CL2D or ML2D.')
        form.addParam('inputClasses', PointerParam, label='Input classes', important=True,
                      pointerClass='SetOfClasses', 
                      help='Select the input classes for the basis construction.')
        
        form.addSection(label='Basis construction')
        form.addParam('maxClasses', IntParam, default=128,
                      label='Max. number of classes', expertLevel=LEVEL_EXPERT,
                      help='Maximum number of classes.')
        form.addParam('maxPCABases', IntParam, default=200,
                      label='Number of PCA bases', expertLevel=LEVEL_EXPERT,
                      help='Number of PCA bases.')
        form.addSection(label='Denoising')
        form.addParam('PCABases2Project', IntParam, default=200,
                      label='Number of PCA bases on which to project', expertLevel=LEVEL_EXPERT,
                      help='Number of PCA bases on which to project.')
        
    def _getDefaultParallel(self):
        """ Return the default value for thread and MPI
        for the parallel section definition.
        """
        return (2, 4)
     
    #--------------------------- INSERT steps functions --------------------------------------------            
    def _insertAllSteps(self):
        """ Insert every step of the protocol"""
        # Create extra dir if needed
        
        
        # Convert input images if necessary
        self._insertFunctionStep('denoiseImages', self.inputParticles.getObjId(), self.inputClasses.getObjId()) 
        
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def denoiseImages(self, inputId, inputClassesId):
        imagesMd = self._getPath('images.xmd')
        writeSetOfParticles(self.inputParticles.get(), imagesMd)
        args = '-i Particles@%s --oroot %s --eigenvectors %d --maxImages %d' % (imagesMd)

        self.runJob("xmipp_image_rotational_pca", args)
        
        args="-i %s -o %s.stk --save_metadata_stack %s.xmd --basis %s.stk %d"
        
        self.runJob("xmipp_transform_filter", args)

        self.outputMd = String(imagesMd)

    def createOutputStep(self):
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(self.inputParticles.get())
        readSetOfParticles(self.outputMd.get(), imgSet)

        self._defineOutputs(outputParticles=imgSet)

    #--------------------------- INFO functions --------------------------------------------                
    def _summary(self):
        import os
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        else:
            summary.append("Fill in with summary information")
        return summary
    
    def _validate(self):
        pass
        
    def _citations(self):
        return []
    
    def _methods(self):
        methods = []
        methods.append('Fill in with methods information.')
        return methods
    
