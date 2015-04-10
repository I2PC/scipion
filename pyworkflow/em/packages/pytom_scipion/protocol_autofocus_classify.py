# **************************************************************************
# *
# * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
# *              J.M. de la Rosa Trevin  (jmdelarosa@cnb.csic.es)
# *              Yuxiang Chen
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

import os

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
import pyworkflow.em as em  

from convert import printLicense, writeSetOfVolumes



class ProtAutofocusClassify(em.ProtClassify3D):
    """ Subtomogram averaging using pytom autofocus """
    _label = 'autofocus'

    #--------------------------- DEFINE param functions --------------------------------------------
    
    def _defineParams(self, form):
        
        form.addSection('Input')
        form.addParam('inputVolumes', params.PointerParam, 
                      pointerClass='SetOfVolumes',
                      label="Input volume particles", important=True, 
                      help='Subtomograms to average')
        form.addParam('provideReferences', params.BooleanParam, default=False,
                      label='Provide initial references?',
                      help='Use _Yes_ if you want to provide you own\n'
                           'set of volumes as initial references.')
        form.addParam('inputReferences', params.PointerParam,
                      condition='provideReferences', allowsNull=True,
                      pointerClass='SetOfVolumes',
                      label='Input reference volumes',
                      help='Set of volumes used as initial references for classification')
        form.addParam('numberOfReferences', params.IntParam, default=4,
                      condition='not provideReferences',
                      label='Number of references', 
                      help="How many references are computed at the end of the process")
        form.addParam('numberOfIterations', params.IntParam, default=10,
                      label='Number of iterations',
                      help='')
        form.addParam('binning', params.IntParam, default=1,
                      label='Binning factor (int)')
        
        form.addSection('Alignment')
        form.addParam('doAlignment', params.BooleanParam, default=True,
                       label='Perform alignment?',
                       help='Set to _No_ for pre-alignned particles.')
        form.addParam('alignmentMask', params.PointerParam,
                       condition='doAlignment',
                       pointerClass='VolumeMask', allowsNull=True,
                       label='Alignment mask',
                       help='Mask used during alignment')
        form.addParam('offset', params.IntParam, 
                      condition='doAlignment', allowsNull=True,
                      label='Alignment offset',
                      help='Maximum number of pixels to shift in alignment.')
        
        form.addSection('Classification')
        form.addParam('maxFreq', params.IntParam, 
                      label='Max. frequency (px)',
                      help='Maximal frequency (in pixels) involved in score calculation.')
        form.addParam('clasificationMask', params.PointerParam,
                       pointerClass='VolumeMask', allowsNull=True,
                       label='Focused classification mask',
                       help='Mask used during classification')
        form.addParam('dispersion', params.IntParam, allowsNull=True,
                       label='Max. cluster dispersion',
                       help='Maximal allowed cluster dispersion.\n'
                            'For example, if you provide a value of 10,\n'
                            'then Nmax / Nmin will always be less than 10.\n'
                            'Where Nmax is the size of biggest cluster and\n'
                            'Nmin the size of smallest cluster.')
        noiseRange = params.Range(0., 1., error='Noise value should be between 0 and 1')
        form.addParam('noise', params.FloatParam, default=0.,
                       validators=[noiseRange],
                       label='Noise value (0 < n < 1)',
                       help='Ratio of the particles to be considered noise')
        line = form.addLine('Difference map',
                             help='    *sigma*: Particle density threshold for difference map.\n'
                                  '*threshold*: STD threshold for difference map.')
        line.addParam('sigma', params.FloatParam, 
                      allowsNull=True, 
                      label='sigma')
        line.addParam('threshold', params.FloatParam, 
                      default=0.4,  allowsNull=True, 
                      label='threshold')
        
        
        form.addParallelSection(threads=0, mpi=4)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):
        """ Mainly prepare the command line for calling reconstruct_significant program"""
        self.volXml = self._getExtraPath('volume_particles.xml')
        
        self._insertFunctionStep('convertInputStep')
        self._insertAutofocusStep()
        self._insertFunctionStep('createOutputStep') 
        
    def _insertAutofocusStep(self):
        """ Prepare the command line to call autofocus_classify.py program. """
        args = '-p %s ' % os.path.basename(self.volXml)
        
        
        if self.provideReferences:
            refList = ','.join('reference_%02d.mrc' % (r+1) for r, _ in enumerate(self.inputReferences.get()))
            args += '-r %s ' % refList
        else:
            args += '-k %d ' % self.numberOfReferences
            
        args += '-i %d ' % self.numberOfIterations
        args += '-b %d ' % self.binning
        
        if self.doAlignment:
            if self.alignmentMask.get() is not None:
                pass #TODO pass mask volume, convert it if needed
            args += '-s %d ' % self.offset if self.offset.hasValue() else ''
        else:
            args += '-a '
            
        args += '-f %d ' % self.maxFreq 
        args += '-d %d ' % self.dispersion if self.dispersion.hasValue() else ''
        args += '-n %0.3f ' % self.noise
        args += '-g %0.3f ' % self.sigma if self.sigma.hasValue() else ''
        args += '-t %0.3f ' % self.threshold if self.threshold.hasValue() else ''
        
        self._insertFunctionStep('runAutofocusStep', args)
        

    #--------------------------- STEPS functions --------------------------------------------        
    def convertInputStep(self):
        printLicense()
        self.info('Wrinting pytom xml file: ' + self.volXml)
        volsDir = self._getExtraPath('inputVolumes')
        pwutils.makePath(volsDir)
        writeSetOfVolumes(self.inputVolumes.get(), self.volXml, volsDir)
        if self.provideReferences:
            ih = em.ImageHandler()
            self.info('Converting input references')
            for r, vol in enumerate(self.inputReferences.get()):
                refFn = self._getExtraPath('reference_%02d.mrc' % (r+1))
                ih.convert(vol, refFn)
            
    def runAutofocusStep(self, args):
        script = self._getScript("classification", "auto_focus_classify.py")
        self.runJob('python', '%s %s' % (script, args), 
                    cwd=self._getExtraPath())

    def createOutputStep(self):
        pass

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        """ Check that some preconditions are met before launching 
        the auto-focus classification run. 
        """
        errors = []
        
        if self.numberOfMpi < 2:
            errors.append('Number of MPI should be greater than 2.')
            
        if self.provideReferences:
            if self.inputReferences.get() is None:
                errors.append('Please select the reference or select *No* in Provide references option')
            
        inputVols = self.inputVolumes.get()
        if inputVols is not None:
            xdim = inputVols.getDim()[0]
            half = xdim/2
            
            if self.maxFreq >= half:
                errors.append('Frequency should be less dim/2 pixels (%d/2=%d)' % (xdim, half))
        
        return errors
        
    def _summary(self):
        summary = []
        return summary
    
    def _citations(self):
        return ['Chen2014']
    
    def _methods(self):
        return []
    
    #--------------------------- UTILS functions --------------------------------------------   
    
    def _getScript(self, *paths):
        return os.path.join(os.environ['PYTOM_HOME'], *paths)
    
    def getDiffMap(self, it, class1, class2):
        return self._getExtraPath('iter%d_dmap_%d_%d.em' % (it, class1, class2))
    
    def getClassMap(self, it, classNo):
        return self._getExtraPath('iter%d_class%d.em' % (it, classNo))
