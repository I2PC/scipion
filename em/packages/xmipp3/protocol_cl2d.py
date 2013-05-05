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
This sub-package contains wrapper around CL2D Xmipp program
"""


from pyworkflow.em import *  
from pyworkflow.utils import *  
import xmipp
from data import *
from xmipp3 import XmippProtocol

# Comparison methods enum
CMP_CORRELATION = 0
CMP_CORRENTROPY = 1

# Clustering methods enum
CL_CLASSICAL = 0
CL_ROBUST = 1


class XmippDefCL2D(Form):
    """Create the definition of parameters for
    the XmippProtCL2D protocol"""
    def __init__(self):
        Form.__init__(self)
    
        self.addSection(label='Input')
        self.addParam('inputImages', PointerParam, label="Input images", important=True, 
                      pointerClass='SetOfImages',
                      help='Select the input images from the project.'
                           'It should be a SetOfImages class')        
        self.addParam('numberOfReferences', IntParam, default=64,
                      label='Number of references:',
                      help='Number of references (or classes) to be generated.')
        self.addParam('numberOfInitialReferences', IntParam, default=4, expertLevel=LEVEL_ADVANCED,
                      label='Number of initial references:',
                      help='Initial number of references used in the first level.')
        self.addParam('numberOfIteration', IntParam, default=4, expertLevel=LEVEL_ADVANCED,
                      label='Number of iterations:',
                      help='Maximum number of iterations within each level.')         
        self.addParam('comparisonMethod', EnumParam, choices=['correlation', 'correntropy'],
                      label="Comparison method", default=CMP_CORRELATION,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Use correlation or correntropy')
        self.addParam('clusteringMethod', EnumParam, choices=['classical', 'robust'],
                      label="Clustering method", default=CL_CLASSICAL,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Use the classical clustering criterion or the robust')
        self.addParam('additionalParams', StringParam, expertLevel=LEVEL_EXPERT,
              label='Additional parameters',
              help='Additional parameters for classify_CL2D:\n  --verbose, --corrSplit, ...')   
        
        self.addSection(label='Core analysis')        
        self.addParam('thZscore', FloatParam, default=3,
                      label='Junk Zscore',
                      help='Which is the average Z-score to be considered as junk. Typical values'
                           'go from 1.5 to 3. For the Gaussian distribution 99.5% of the data is'
                           'within a Z-score of 3. Lower Z-scores reject more images. Higher Z-scores'
                           'accept more images.')
        self.addParam('thPCAZscore', FloatParam, default=3,
                      label='PCA Zscore',
                      help='Which is the PCA Z-score to be considered as junk. Typical values'
                           'go from 1.5 to 3. For the Gaussian distribution 99.5% of the data is'
                           'within a Z-score of 3. Lower Z-scores reject more images. Higher Z-scores'
                           'accept more images.')        
        self.addParam('tolerance', IntParam, default=1,
                      label='Tolerance',
                      help='An image belongs to the stable core if it has been with other images in the same class'
                           'in all the previous levels except possibly a few of them. Tolerance defines how few is few.'
                           'Tolerance=0 means that an image must be in all previous levels with the rest of images in'
                           'the core.')          
        
        self.addParallelSection(threads=0, mpi=2)
        
        
class XmippProtCL2D(ProtAlign, ProtClassify, XmippProtocol):
    """ Protocol to preprocess a set of micrographs in the project. """
    _definition = XmippDefCL2D()
    _label = 'Xmipp CL2D'

    def _defineSteps(self):
        """ Mainly prepare the command line for call ml(f)2d program"""
        progId = self.getProgramId()
        
        program = "xmipp_%s_align2d" % progId

        restart = False
        
        prefix = '%s2d_' % progId
        self.oroot = self._getPath(prefix)
        
        self.inputImgs = self.inputImages.get()
        
        imgsFn = self._insertConvertStep('inputImgs', XmippSetOfImages,
                                         self._getPath('input_images.xmd'))
        params = ' -i %s --oroot %s' % (imgsFn, self.oroot)
        # Number of references will be ignored if -ref is passed as expert option
        if self.doGenerateReferences:
            params += ' --nref %d' % self.numberOfReferences.get()
        else:
            self.inputRefs = self.inputReferences.get()
            refsFn = self._insertConvertStep('inputRefs', XmippSetOfImages,
                                             self._getPath('input_references.xmd'))
            params += ' --ref %s' % refsFn
        
        if self.doMlf:
            if not self.doCorrectAmplitudes:
                params += ' --no_ctf'                    
            if not self.areImagesPhaseFlipped:
                params += ' --not_phase_flipped'
            if self.highResLimit.get() > 0:
                params += ' --limit_resolution 0 %f' % self.highResLimit.get()
            params += ' --sampling_rate %f' % self.inputImages.get().samplingRate.get()
        else:
            if self.doFast:
                params += ' --fast'
            if self.numberOfThreads.get() > 1:
                params += ' --thr %d' % self.numberOfThreads.get()
            
        if self.maxIters.get() != 100:
            params += " --iter %d" % self.maxIters.get()

        if self.doMirror:
            params += ' --mirror'
            
        if self.doNorm:
            params += ' --norm'

        self._insertRunJobStep(program, params)
                
        self._insertFunctionStep('createOutput')
        
    def createOutput(self):
        classification = XmippClassification2D(self.oroot + 'classes.xmd')
        self._defineOutputs(outputClassification=classification)
