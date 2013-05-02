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
This sub-package contains wrapper around ML2D Xmipp program
"""


from pyworkflow.em import *  
from pyworkflow.utils import *  
import xmipp
from data import *
from xmipp3 import XmippProtocol


class XmippDefML2D(Form):
    """Create the definition of parameters for
    the XmippProtML2D protocol"""
    def __init__(self):
        Form.__init__(self)
    
        self.addSection(label='Input')
        self.addParam('inputImages', PointerParam, label="Input images", important=True, 
                      pointerClass='SetOfImages',
                      help='Select the input images from the project.'
                           'It should be a SetOfImages class')        
        self.addParam('doGenerateReferences', BooleanParam, default=True,
                      label='Generate references (or classes) ?', 
                      help='If you set to <No>, you should provide references images'
                           'If <Yes>, the default generation is done by averaging'
                           'subsets of the input images.')
        self.addParam('numberOfReferences', IntParam, default=3, condition='doGenerateReferences',
                      label='Number of references:',
                      help='Number of references to be generated.')
        self.addParam('referenceImages', PointerParam, condition='not doGenerateReferences',
                      label="Reference image(s)", 
                      pointerClass='SetOfImages',
                      help='Image(s) that will serve as class references')
        
        self.addSection(label='MLF-specific parameters', questionParam='doMlf')        
        self.addParam('doMlf', BooleanParam, default=False,
                      label='Use MLF2D instead of ML2D')
        self.addParam('doCorrectAmplitudes', BooleanParam, default=True,
                      label='Use CTF-amplitude correction inside MLF?',
                      help='If set to <Yes>, the input images file should contains'
                           'the CTF information for each image.'
                           'If set to <No>, provide the images pixel size in Angstrom.')
        self.addParam('areImagesPhaseFlipped', BooleanParam, default=True,
                      label='Are the images CTF phase flipped?',
                      help='You can run MLF with or without having phase flipped the images.')        
        self.addParam('highResLimit', IntParam, default=20,
                      label='High-resolution limit (Angstroms)',
                      help='No frequencies higher than this limit will be taken into account.'
                           'If zero is given, no limit is imposed.')
        
        self.addSection(label='Advanced parameters', questionParam='showAdvanced')        
        self.addParam('showAdvanced', BooleanParam, default=False,
                      label='Show advanced parameters')
        self.addParam('doMirror', BooleanParam, default=True,
                      label='Also include mirror in the alignment?',
                      help='Including the mirror transformation is useful if your particles'
                           'have a handedness and may fall either face-up or face-down on the grid.'
                           )
        self.addParam('doFast', BooleanParam, default=True, condition='not doMlf',
                      label='Use the fast version of this algorithm?',
                      help='For details see:\n'
                           '<Scheres et al., Bioinformatics, 21 (Suppl. 2), ii243-ii244>\n'
                           '[http://dx.doi.org/10.1093/bioinformatics/bti1140]'
                           )        
        self.addParam('doNorm', BooleanParam, default=False,
                      label='Refine the normalization for each image?',
                      help='This variant of the algorithm deals with normalization errors.'
                           'For more info see (and please cite):\n'
                           '<Scheres et. al. (2009) J. Struc. Biol., Vol 166, Issue 2, May 2009>\n'
                           '[http://dx.doi.org/10.1016/j.jsb.2009.02.007]'
                           )             
        # Advance or expert parameters
        self.addParam('maxIters', IntParam, default=100, expertLevel=LEVEL_ADVANCED,
                      label='Maximum number of iterations',
                      help='If the convergence has not been reached after this number'
                           'of iterations, the process will be stopped.')   
        self.addParam('psiStep', FloatParam, default=5.0, expertLevel=LEVEL_ADVANCED,
                      label='In-plane rotation sampling (degrees)',
                      help='In-plane rotation sampling interval (degrees).')          
        self.addParam('stdNoise', FloatParam, default=1.0, expertLevel=LEVEL_EXPERT,
                      label='Std for pixel noise',
                      help='Expected standard deviation for pixel noise.')               
        self.addParam('stdOffset', FloatParam, default=3.0, expertLevel=LEVEL_EXPERT,
                      label='Std for origin offset',
                      help='Expected standard deviation for origin offset (pixels).') 
        
        self.addParam('numberOfThreads', IntParam, default=2, 
                      label='Maximum number of threads',
                      help='If the convergence has not been reached after this number'
                           'of iterations, the process will be stopped.')          
        
        
class XmippProtML2D(ProtAlign, ProtClassify, XmippProtocol):
    """ Protocol to preprocess a set of micrographs in the project. """
    _definition = XmippDefML2D()
    _label = 'Xmipp ML2D'

    def getProgramId(self):
        progId = "ml"
        if self.doMlf:
            progId += "f" 
        return progId   
         
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
