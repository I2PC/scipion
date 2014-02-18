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

        
class XmippProtML2D(ProtClassify):
    """
    Perform (multi-reference) 2D-alignment using 
    a maximum-likelihood ( *ML* ) target function.
    
    Initial refereces can be generated from random subsets of the experimental 
    images or can be provided by the user (this can introduce bias). The output 
    of the protocol consists of the refined 2D classes (weighted averages over 
    all experimental images). The experimental images are not altered at all.    
    
    Although the calculations can be rather time-consuming (especially for 
    many, large experimental images and a large number of references we 
    strongly recommend to let the calculations converge. 
    """
    _label = 'ml2d'
    
    def __init__(self, **args):
        ProtClassify.__init__(self, **args)        
        self.progId = "ml"
        self.oroot = ""
       
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam, label="Input particles", important=True, 
                      pointerClass='SetOfParticles',
                      help='Select the input images from the project.')        
        form.addParam('doGenerateReferences', BooleanParam, default=True,
                      label='Generate references?', 
                      help='If you set to *No*, you should provide references images'
                           'If *Yes*, the default generation is done by averaging'
                           'subsets of the input images. (less bias introduced)')
        form.addParam('numberOfReferences', IntParam, default=3, condition='doGenerateReferences',
                      label='Number of references:',
                      help='Number of references to be generated.')
        form.addParam('referenceParticles', PointerParam, condition='not doGenerateReferences',
                      label="Reference image(s)", 
                      pointerClass='SetOfParticles',
                      help='Image(s) that will serve as initial 2D references')
        
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
                      label='Use the fast version?',
                      help='If set to *Yes*, a fast approach will be used to avoid\n'
                           'searching in the whole solutions space.             \n\n'
                           'For details see (and please cite): \n' + self._getCite('Scheres2005b')
                           )
        form.addParam('doNorm', BooleanParam, default=False,
                      label='Refine the normalization for each image?',
                      help='This variant of the algorithm deals with normalization errors. \n\n'
                           'For details see (and please cite): \n ' + self._getCite('Scheres2009b')
                           )             
        # Advance or expert parameters
        form.addParam('maxIters', IntParam, default=100, expertLevel=LEVEL_ADVANCED,
                      label='Maximum number of iterations',
                      help='If the convergence has not been reached after this number'
                           'of iterations, the process will be stopped.')   
        form.addParam('psiStep', FloatParam, default=5.0, expertLevel=LEVEL_ADVANCED,
                      label='In-plane rotation sampling (degrees)',
                      help='In-plane rotation sampling interval (degrees).')          
        form.addParam('stdNoise', FloatParam, default=1.0, expertLevel=LEVEL_EXPERT,
                      label='Std for pixel noise',
                      help='Expected standard deviation for pixel noise.')               
        form.addParam('stdOffset', FloatParam, default=3.0, expertLevel=LEVEL_EXPERT,
                      label='Std for origin offset',
                      help='Expected standard deviation for origin offset (pixels).') 
        
        form.addParallelSection(threads=2, mpi=2)
           
    #--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self):        
        self.oroot = self._getOroot()
        self.program = "xmipp_%s_align2d" % self.progId       
        
        # Convert input images if necessary
        self.imgsFn = createXmippInputImages(self, self.inputParticles.get())
        
        self._insertMLStep()
        self._insertFunctionStep('createOutputStep')

    def _insertMLStep(self):
        """ Mainly prepare the command line for call ml(f)2d program"""
        params = ' -i %s --oroot %s' % (self.imgsFn, self.oroot)
        if self.doGenerateReferences:
            params += ' --nref %d' % self.numberOfReferences.get()
        else:
            self.inputRefs = self.referenceParticles.get()
            refsFn = createXmippInputImages(self, self.inputRefs, imagesFn='input_references.xmd')
            params += ' --ref %s' % refsFn
            self.numberOfReferences.set(self.inputRefs.getSize())
        
        if self.doMlf:
            if not self.doCorrectAmplitudes:
                params += ' --no_ctf'                    
            if not self.areImagesPhaseFlipped:
                params += ' --not_phase_flipped'
            if self.highResLimit > 0:
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
        
    #--------------------------- STEPS functions --------------------------------------------       
    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        classes2DSet = self._createSetOfClasses2D(imgSet)
        readSetOfClasses2D(classes2DSet, self.oroot + 'classes.xmd')
        self._defineOutputs(outputClasses=classes2DSet)
        self._defineSourceRelation(imgSet, classes2DSet)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        cites = ['Scheres2005a']
        
        if self.doMlf:
            cites.append('Scheres2007b')
        elif self.doFast:
            cites.append('Scheres2005b')
            
        if self.doNorm:
            cites.append('Scheres2009b')
            
        return cites
    
    def _summary(self):
        summary = []
        summary.append('Number of input images: *%d*' % self.inputParticles.get().getSize())
        summary.append('Classified into *%d* classes' % self.numberOfReferences.get())
        
        if self.doMlf:
            summary.append('- Used a ML in _Fourier-space_')
        elif self.doFast:
            summary.append('- Used _fast_, reduced search-space approach')

        if self.doNorm:
            summary.append('- Refined _normalization_ for each experimental image')
            
        return summary
    
    def _methods(self):
        return self._summary()  # summary is quite explicit and serve as methods
    
    #--------------------------- UTILS functions --------------------------------------------
    def _getIterClasses(self, it=None, block=None):
        """ Return the classes metadata for this iteration.
        block parameter can be 'info' or 'classes'.
        """
        if it is None:
            it = self._lastIteration()
        extra = self.oroot + 'extra'
        mdFile = join(extra, 'iter%03d' % it, 'iter_classes.xmd')
        if block:
            mdFile = block + '@' + mdFile
        
        return mdFile
    
    def _lastIteration(self):
        """ Find the last iteration number """
        if self.oroot == "":
            self.oroot = self._getOroot()
        it = 0        
        while True:
            if not exists(self._getIterClasses(it+1)):
                break
            it += 1
        return it        
    
    def _getOroot(self):        
        if self.doMlf:
            self.progId += "f"
        return self._getPath('%s2d_' % self.progId)  
    