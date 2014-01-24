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
from convert import createXmippInputImages, readSetOfClasses2D, readSetOfParticles

#from xmipp3 import XmippProtocol
        
        
class XmippProtML2D(ProtClassify):
    """ Protocol to preprocess a set of micrographs in the project. """
    _label = 'ml2d'
    
    def __init__(self, **args):
        ProtClassify.__init__(self, **args)
        
        self.progId = "ml"
        self.oroot = ""
        
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam, label="Input particles", important=True, 
                      pointerClass='SetOfParticles',
                      help='Select the input images from the project.'
                           'It should be a SetOfImages class')        
        form.addParam('doGenerateReferences', BooleanParam, default=True,
                      label='Generate references?', 
                      help='If you set to *No*, you should provide references images'
                           'If *Yes*, the default generation is done by averaging'
                           'subsets of the input images.')
        form.addParam('numberOfReferences', IntParam, default=3, condition='doGenerateReferences',
                      label='Number of references:',
                      help='Number of references to be generated.')
        form.addParam('referenceParticles', PointerParam, condition='not doGenerateReferences',
                      label="Reference image(s)", 
                      pointerClass='SetOfImages',
                      help='Image(s) that will serve as class references')
        
        #form.addSection(label='MLF-specific parameters', questionParam='doMlf')        
        form.addParam('doMlf', BooleanParam, default=False, important=True,
                      label='Use MLF2D instead of ML2D?')
        form.addParam('doCorrectAmplitudes', BooleanParam, default=True, condition='doMlf',
                      label='Use CTF-amplitude correction?',
                      help='If set to *Yes*, the input images file should contains'
                           'If set to *No*, provide the images pixel size in Angstrom.')
        form.addParam('areImagesPhaseFlipped', BooleanParam, default=True, condition='doMlf',
                      label='Are the images CTF phase flipped?',
                      help='You can run MLF with or without having phase flipped the images.')        
        form.addParam('highResLimit', IntParam, default=20, condition='doMlf',
                      label='High-resolution limit (Ang)',
                      help='No frequencies higher than this limit will be taken into account.'
                           'If zero is given, no limit is imposed.')
        
        form.addSection(label='Advanced')#, questionParam='showAdvanced')        
#        form.addParam('showAdvanced', BooleanParam, default=False,
#                      label='Show advanced parameters')
        form.addParam('doMirror', BooleanParam, default=True,
                      label='Also include mirror in the alignment?',
                      help='Including the mirror transformation is useful if your particles'
                           'have a handedness and may fall either face-up or face-down on the grid.'
                           )
        form.addParam('doFast', BooleanParam, default=True, condition='not doMlf',
                      label='Use the fast version of this algorithm?',
                      help='For details see (and please cite): \n '
                           '*Scheres et al., Bioinformatics, 21 (Suppl. 2), ii243-ii244* \n '
                           '[[http://dx.doi.org/10.1093/bioinformatics/bti1140][Info]]'
                           )        
        form.addParam('doNorm', BooleanParam, default=False,
                      label='Refine the normalization for each image?',
                      help='This variant of the algorithm deals with normalization errors. \n '
                           'For details see (and please cite): \n '
                           '*Scheres et. al. (2009) J. Struc. Biol., Vol 166, Issue 2, May 2009* \n '
                           '[[http://dx.doi.org/10.1016/j.jsb.2009.02.007][Info]]'
                           )             
        # Advance or expert parameters
        form.addParam('maxIters', IntParam, default=100,# expertLevel=LEVEL_ADVANCED,
                      label='Maximum number of iterations',
                      help='If the convergence has not been reached after this number'
                           'of iterations, the process will be stopped.')   
        form.addParam('psiStep', FloatParam, default=5.0,# expertLevel=LEVEL_ADVANCED,
                      label='In-plane rotation sampling (degrees)',
                      help='In-plane rotation sampling interval (degrees).')          
        form.addParam('stdNoise', FloatParam, default=1.0, expertLevel=LEVEL_EXPERT,
                      label='Std for pixel noise',
                      help='Expected standard deviation for pixel noise.')               
        form.addParam('stdOffset', FloatParam, default=3.0, expertLevel=LEVEL_EXPERT,
                      label='Std for origin offset',
                      help='Expected standard deviation for origin offset (pixels).') 
        
        form.addParallelSection(threads=2, mpi=2)
            
    def _getIterClasses(self, iter=None, block=None):
        """ Return the classes metadata for this iteration.
        block parameter can be 'info' or 'classes'.
        """
        if iter is None:
            iter = self._lastIteration()
        extra = self.oroot + 'extra'
        mdFile = join(extra, 'iter%03d' % iter, 'iter_classes.xmd')
        if block:
            mdFile = block + '@' + mdFile
        
        return mdFile
    
    def _lastIteration(self):
        """ Find the last iteration number """
        if self.oroot == "":
            self.oroot = self._getOroot()
        iterNumber = 0        
        while True:
            if not exists(self._getIterClasses(iterNumber+1)):
                break
            iterNumber = iterNumber + 1
        return iterNumber        
    
    def _getOroot(self):
        
        if self.doMlf:
            self.progId += "f"
        return self._getPath('%s2d_' % self.progId)       
        
    def _defineSteps(self):
        """ Mainly prepare the command line for call ml(f)2d program"""
        
        self.oroot = self._getOroot()
        self.program = "xmipp_%s_align2d" % self.progId       
        
        # Convert input images if necessary
        imgsFn = createXmippInputImages(self, self.inputParticles.get())
        
        params = ' -i %s --oroot %s' % (imgsFn, self.oroot)
        # Number of references will be ignored if -ref is passed as expert option
        if self.doGenerateReferences:
            params += ' --nref %d' % self.numberOfReferences.get()
        else:
            self.inputRefs = self.referenceParticles.get()
            refsFn = createXmippInputImages(self, self.inputRefs, imagesFn='input_references.xmd')
            params += ' --ref %s' % refsFn
        
        if self.doMlf:
            if not self.doCorrectAmplitudes:
                params += ' --no_ctf'                    
            if not self.areImagesPhaseFlipped:
                params += ' --not_phase_flipped'
            if self.highResLimit.get() > 0:
                params += ' --limit_resolution 0 %f' % self.highResLimit.get()
            params += ' --sampling_rate %f' % self.inputParticles.get().getSamplingRate()
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

        self._insertRunJobStep(self.program, params)
                
        self._insertFunctionStep('createOutput')
        
    def createOutput(self):
        
        classes2DSet = self._createSetOfClasses2D()
        classes2DSet.setImages(self.inputParticles.get())
        readSetOfClasses2D(classes2DSet, self.oroot + 'classes.xmd')
        print "after read"
        self._defineOutputs(outputClasses=classes2DSet)
        print "after define"
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputClasses'):
            summary.append("Output classes not ready yet.")
        else:
            summary.append("Input Images: %s" % self.inputParticles.get().getNameId())
            summary.append("Number of references: %d" % self.numberOfReferences.get())
            summary.append("Output classes: %s" % self.outputClasses.get())
        return summary