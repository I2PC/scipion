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

from os.path import join, exists

from pyworkflow.em.protocol import ProtClassify2D
from pyworkflow.protocol.constants import LEVEL_ADVANCED, LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, BooleanParam, IntParam, FloatParam

from convert import writeSetOfParticles, readSetOfClasses2D


        
class XmippProtML2D(ProtClassify2D):
    """
    Perform (multi-reference) 2D-alignment using 
    a maximum-likelihood ( *ML* ) target function.
    
    Initial references can be generated from random subsets of the experimental 
    images or can be provided by the user (this can introduce bias). The output 
    of the protocol consists of the refined 2D classes (weighted averages over 
    all experimental images). The experimental images are not altered at all.    
    
    Although the calculations can be rather time-consuming (especially for 
    many, large experimental images and a large number of references we 
    strongly recommend to let the calculations converge. 
    """
    _label = 'ml2d'
    
    def __init__(self, **args):
        ProtClassify2D.__init__(self, **args)
        
    def _defineFileNames(self):
        """ Centralize how files are called within the protocol. """
        myDict = {
                  'input_particles': self._getTmpPath('input_particles.xmd'),
                  'input_references': self._getTmpPath('input_references.xmd'),
                  'output_classes': self._getOroot() + 'classes.xmd',
                  }
        self._updateFilenamesDict(myDict)
       
    #--------------------------- DEFINE param functions --------------------------------------------   
    
    def _defineParams(self, form):
        form.addSection(label='Params')
        
        group = form.addGroup('Input')
        group.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles',
                       label="Input particles", important=True,
                       help='Select the input images from the project.')        
        group.addParam('doGenerateReferences', BooleanParam, default=True,
                      label='Generate classes?',
                      help='If you set to *No*, you should provide class images.\n'
                           'If *Yes*, the default generation is done by averaging\n'
                           'subsets of the input images (less bias introduced).')
        group.addParam('numberOfClasses', IntParam, default=3, condition='doGenerateReferences',
                      label='Number of classes:',
                      help='Number of classes to be generated.')
        group.addParam('inputReferences', PointerParam, allowsNull=True,
                       condition='not doGenerateReferences',
                      label="Class image(s)",
                      pointerClass='SetOfParticles',
                      help='Image(s) that will serve as initial 2D classes')
        
        form.addParam('doMlf', BooleanParam, default=False, important=True,
                      label='Use MLF2D instead of ML2D?')
        
        group = form.addGroup('ML-Fourier', condition='doMlf')
        group.addParam('doCorrectAmplitudes', BooleanParam, default=True,
                      label='Use CTF-amplitude correction?',
                      help='If set to *Yes*, the input images file should contains.\n'
                           'If set to *No*, provide the images pixel size in Angstrom.')
        group.addParam('areImagesPhaseFlipped', BooleanParam, default=True,
                      label='Are the images CTF phase flipped?',
                      help='You can run MLF with or without having phase flipped the images.')        
        group.addParam('highResLimit', IntParam, default=20,
                      label='High-resolution limit (Ang)',
                      help='No frequencies higher than this limit will be taken into account.\n'
                           'If zero is given, no limit is imposed.')
        
        form.addSection(label='Advanced')
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
        form.addParam('stdNoise', FloatParam, default=1.0, expertLevel=LEVEL_ADVANCED,
                      label='Std for pixel noise',
                      help='Expected standard deviation for pixel noise.')               
        form.addParam('stdOffset', FloatParam, default=3.0, expertLevel=LEVEL_ADVANCED,
                      label='Std for origin offset',
                      help='Expected standard deviation for origin offset (pixels).') 
        
        form.addParallelSection(threads=2, mpi=4)
           
    #--------------------------- INSERT steps functions --------------------------------------------  
    
    def _insertAllSteps(self):  
        self._defineFileNames()      
        self._insertFunctionStep('convertInputStep', self.inputParticles.get().getObjId())
        program = self._getMLProgram()
        params = self._getMLParams()
        self._insertRunJobStep(program, params)
        self._insertFunctionStep('createOutputStep')

    def _getMLParams(self):
        """ Mainly prepare the command line for call ml(f)2d program"""
        params = ' -i %s --oroot %s' % (self._getFileName('input_particles'), self._getOroot())
        if self.doGenerateReferences:
            params += ' --nref %d' % self.numberOfClasses.get()
            self.inputReferences.set(None)
        else:
            params += ' --ref %s' % self._getFileName('input_references')
            self.numberOfClasses.set(self.inputReferences.get().getSize())
        
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
            if self.numberOfThreads > 1:
                params += ' --thr %d' % self.numberOfThreads.get()
            
        if self.maxIters != 100:
            params += ' --iter %d' % self.maxIters.get()

        if self.doMirror:
            params += ' --mirror'
            
        if self.doNorm:
            params += ' --norm'

        return params
        
    #--------------------------- STEPS functions --------------------------------------------       
    def convertInputStep(self, inputId):
        """ Write the input images as a Xmipp metadata file. """
        writeSetOfParticles(self.inputParticles.get(), self._getFileName('input_particles'))
        # If input references, also convert to xmipp metadata
        if not self.doGenerateReferences:
            writeSetOfParticles(self.inputReferences.get(), self._getFileName('input_references'))
        
    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        classes2DSet = self._createSetOfClasses2D(imgSet)
        readSetOfClasses2D(classes2DSet, self._getFileName('output_classes'))
        self._defineOutputs(outputClasses=classes2DSet)
        self._defineSourceRelation(self.inputParticles, classes2DSet)
        if not self.doGenerateReferences:
            self._defineSourceRelation(self.inputReferences, classes2DSet)
    
    #--------------------------- INFO functions -------------------------------------------- 
    
    def _validate(self):
        errors = []
        if self.doMlf:
            inputParticles = self.inputParticles.get()
            if inputParticles is not None and not inputParticles.hasCTF():
                errors.append('Input particles does not have CTF information.\n'
                              'This is required when using ML in fourier space.')
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
        if hasattr(self, 'outputClasses'):
            summary.append('Input Particles: *%d*' % self.inputParticles.get().getSize())
            summary.append('Classified into *%d* classes' % self.numberOfClasses.get())

            if self.doMlf:
                summary.append('- Used a ML in _Fourier-space_')
            elif self.doFast:
                summary.append('- Used _fast_, reduced search-space approach')

            if self.doNorm:
                summary.append('- Refined _normalization_ for each experimental image')
            
        return summary
    
    def _methods(self):
        methods = []
        if hasattr(self, 'outputClasses'):
            methods.append('Input dataset %s of *%d* images was classified' % (self.getObjectTag('inputParticles'), self.inputParticles.get().getSize()))
            numberOfClasses = self.numberOfClasses.get()
            classesTxt =  'class' if numberOfClasses == 1 else 'classes'
            methods.append('into *%d* 2D %s using Maximum Likelihood (ML) inside Xmipp.' % (numberOfClasses, classesTxt))

            if self.doMlf:
                methods.append('ML was used in _Fourier-space_.')
            elif self.doFast:
                methods.append('Used _fast_, reduced search-space approach.')

            if self.doNorm:
                methods.append('The _normalization_ was refined for each experimental image.')
            methods.append('Output set is %s.'%(self.getObjectTag('outputClasses')))

        return methods

    #--------------------------- UTILS functions --------------------------------------------
    
    def _getIterClasses(self, it=None, block=None):
        """ Return the classes metadata for this iteration.
        block parameter can be 'info' or 'classes'.
        """
        if it is None:
            it = self._lastIteration()
            
        extra = self._getOroot() + 'extra'
        mdFile = join(extra, 'iter%03d' % it, 'iter_classes.xmd')
        if block:
            mdFile = block + '@' + mdFile
        
        return mdFile
    
    def _lastIteration(self):
        """ Find the last iteration number """
        it = 0        
        while True:
            if not exists(self._getIterClasses(it+1)):
                break
            it += 1
        return it        
    
    def _getMLId(self):
        """ Return ml or mlf depending if using fourier or not. """
        if self.doMlf:
            return 'mlf'
        return 'ml'
    
    def _getMLProgram(self):
        """ Return the program to be used, depending if using fourier. """
        return "xmipp_%s_align2d" % self._getMLId()
    
    def _getOroot(self):        
        return self._getPath('%s2d_' % self._getMLId())  
    