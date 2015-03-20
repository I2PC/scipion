# **************************************************************************
# *
# * Authors:     Airen Zaldivar Peraza (azaldivar@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
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
In this module are protocol base classes related to EM Micrographs
"""

import os
from os.path import join, basename, exists, dirname, relpath
from itertools import izip

from pyworkflow.object import String, Boolean
from pyworkflow.protocol.constants import STEPS_PARALLEL, LEVEL_ADVANCED, LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, FloatParam, IntParam, TextParam, BooleanParam, FileParam
from pyworkflow.utils.path import copyTree, copyFile, removeBaseExt, makePath, moveFile
from pyworkflow.utils.properties import Message
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.data import Micrograph, SetOfImages, SetOfCTF



class ProtMicrographs(EMProtocol):
    pass


class ProtCTFMicrographs(ProtMicrographs):
    """ Base class for all protocols that estimates the CTF"""
    
    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL
        self.isFirstTime = Boolean(False)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_CTF_ESTI)
        form.addHidden('recalculate', BooleanParam, default=False)
        
        form.addParam('continueRun', PointerParam, allowsNull=True,
                      condition='recalculate', label="input previous run",
                      pointerClass=self.getClassName())
        form.addHidden('sqliteFile', FileParam, condition='recalculate',
                       allowsNull=True)
        
        form.addParam('inputMicrographs', PointerParam, important=True,
                       condition='not recalculate', label=Message.LABEL_INPUT_MIC,
                       pointerClass='SetOfMicrographs')
        form.addParam('ctfDownFactor', FloatParam, default=1.,
                      label='CTF Downsampling factor',
                      condition='not recalculate',
                      help='Set to 1 for no downsampling. Non-integer downsample factors are possible. '
                      'This downsampling is only used for estimating the CTF and it does not affect '
                      'any further calculation. Ideally the estimation of the CTF is optimal when '
                      'the Thon rings are not too concentrated at the origin (too small to be seen) '
                      'and not occupying the whole power spectrum (since this downsampling might '
                      'entail aliasing).')
        
        self._defineProcessParams(form)

        line = form.addLine('Resolution', condition='not recalculate',
                            help='Give a value in digital frequency (i.e. between 0.0 and 0.5). '
                                 'These cut-offs prevent the typical peak at the center of the PSD and high-resolution'
                                 'terms where only noise exists, to interfere with CTF estimation. The default lowest '
                                 'value is 0.05 but for micrographs with a very fine sampling this may be lowered towards 0.'
                                 'The default highest value is 0.35, but it should '+'be increased for micrographs with '
                                 'signals extending beyond this value. However, if your micrographs extend further than '
                                 '0.35, you should consider sampling them at a finer rate.')
        line.addParam('lowRes', FloatParam, default=0.05,
                      label='Lowest' )
        line.addParam('highRes', FloatParam, default=0.35,
                      label='Highest')
        # Switched (microns) by 'in microns' by fail in the identifier with jquery
        line = form.addLine('Defocus search range (microns)', expertLevel=LEVEL_ADVANCED,
                            condition='not recalculate',
                            help='Select _minimum_ and _maximum_ values for defocus search range (in microns).'
                                 'Underfocus is represented by a positive number.')
        line.addParam('minDefocus', FloatParam, default=0.25, 
                      label='Min')
        line.addParam('maxDefocus', FloatParam, default=4.,
                      label='Max')
        
        form.addParam('windowSize', IntParam, default=256, expertLevel=LEVEL_ADVANCED,
                      label='Window size', condition='not recalculate',
                      help='The PSD is estimated from small patches of this size. Bigger patches '
                           'allow identifying more details. However, since there are fewer windows, '
                           'estimations are noisier.')
        
        form.addParallelSection(threads=2, mpi=1)       
    
    def _defineProcessParams(self, form):
        """ This method should be implemented by subclasses
        to add other parameter relatives to the specific operation."""
        pass
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ Insert the steps to perform CTF estimation, or re-estimation, on a set of micrographs.
        """
        deps = [] # Store all steps ids, final step createOutput depends on all of them
        fDeps = []
        
        if not self.recalculate:
            deps = self._insertEstimationSteps()
            # Insert step to create output objects
            fDeps = self._insertFinalSteps(deps)
        else:
            if self.isFirstTime:
                self._insertPreviousSteps() # Insert previous estimation or re-estimation an so on...
                self.isFirstTime.set(False)
            fDeps = self._insertRecalculateSteps()
        
        self._insertFunctionStep('createOutputStep', prerequisites=fDeps)
    
    def _insertFinalSteps(self, deps):
        """ This should be implemented in subclasses"""
        return deps
    
    def _insertEstimationSteps(self):
        estimDeps = []
        self._defineValues()
        self._prepareCommand()
        # For each micrograph insert the steps to process it
        for micFn, micDir, _ in self._iterMicrographs():
            # CTF estimation
            # Make estimation steps independent between them
            stepId = self._insertFunctionStep('_estimateCTF', micFn, micDir,
                                                  prerequisites=[]) # Make estimation steps independent between them
            estimDeps.append(stepId)
        return estimDeps
    
    def _insertRecalculateSteps(self):
        recalDeps = []
        # For each psd insert the steps to process it
        self.recalculateSet = SetOfCTF(filename=self.sqliteFile.get(), objDoStore=False)
        for ctf in self.recalculateSet:
            line = ctf.getObjComment()
            if ctf.isEnabled() and line:
                # CTF Re-estimation
                copyId = self._insertFunctionStep('copyMicDirectoryStep', ctf.getObjId())
                # Make estimation steps independent between them
                stepId = self._insertFunctionStep('_restimateCTF', ctf.getObjId(), prerequisites=[copyId])
                recalDeps.append(stepId)
        return recalDeps
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _estimateCTF(self, micFn, micDir):
        """ Do the CTF estimation with the specific program
        and the parameters required.
        Params:
         micFn: micrograph filename
         micDir: micrograph directory
        """
        raise Exception(Message.ERROR_NO_EST_CTF)
    
    def _restimateCTF(self, id):
        """ Do the CTF estimation with the specific program
        and the parameters required.
        Params:
         micFn: micrograph filename
         micDir: micrograph directory
        """
        raise Exception(Message.ERROR_NO_EST_CTF)
    
    def copyMicDirectoryStep(self, id):
        """ Copy micrograph's directory tree for recalculation"""
        ctfModel = self.recalculateSet[id]
        mic = ctfModel.getMicrograph()
        
        prevDir = self._getPrevMicDir(ctfModel)
        micDir = self._getMicrographDir(mic)
        if not prevDir == micDir:
            # Create micrograph dir under extra directory
            makePath(micDir)
            if not exists(micDir):
                raise Exception("No created dir: %s " % micDir)
            copyTree(prevDir, micDir)
    
    def _createNewCtfModel(self, mic):
        """ This should be implemented in subclasses
        in order to create a CTF model 
        """
        pass
    
    def createOutputStep(self):
        """ This function is shared by Xmipp and CTFfind
        estimation, or recalculate, protocols.
        if is recalculate, it will iterated for each CTF model, see
        if was recalculated and update with new defocus values.
        Else, the function that should be implemented in each subclass.
        """
        if self.recalculate:
            ctfSet = self._createSetOfCTF("_recalculated")
            defocusList = []
            if self.continueRun:
                oldCtfSet = getattr(self.continueRun.get(), 'outputCTF')
            else:
                oldCtfSet = getattr(self, 'outputCTF')
            micSet = oldCtfSet.getMicrographs()
            # README: We suppose this is reading the ctf selection (with enabled/disabled)
            # to only consider the enabled ones in the final SetOfCTF
            
            #TODO: maybe we can remove the need of the extra text file
            # with the recalculate parameters
            for ctfModel in self.recalculateSet:
                if ctfModel.isEnabled() and ctfModel.getObjComment():
                    mic = ctfModel.getMicrograph()
                    # Update the CTF models that where recalculated
                    # and append later to the set
                    # we dont want to copy the id here since it is already correct
                    ctfModel.copy(self._createNewCtfModel(mic), copyId=False)
                    ctfModel.setEnabled(True)
                ctfSet.append(ctfModel)
                # save the values of defocus for each micrograph in a list
                defocusList.append(ctfModel.getDefocusU())
                defocusList.append(ctfModel.getDefocusV())
            ctfSet.setMicrographs(micSet)
            self._defineOutputs(outputCTF=ctfSet)
            self._defineCtfRelation(micSet, ctfSet)
    
            self._defocusMaxMin(defocusList)
            self._ctfCounter(defocusList)
        else:
            self._createOutputStep()
        
    #--------------------------- INFO functions ----------------------------------------------------
    def _summary(self):
        summary = []
        
        if self.recalculate:
            if self.isFinished():
                if self.summaryVar.hasValue():
                    summary.append(self.summaryVar.get())
            else:
                summary.append(Message.TEXT_NO_CTF_READY)
        else:
            if not hasattr(self, 'outputCTF'):
                summary.append(Message.TEXT_NO_CTF_READY)
            else:
                summary.append("CTF estimation of %d micrographs." % self.inputMicrographs.get().getSize())
        
        return summary
    
    def _methods(self):
        methods = []
        
        if hasattr(self, 'outputCTF') and self.isFinished():
            methods.append(self.methodsVar.get())
        else:
            methods.append(Message.TEXT_NO_CTF_READY)
            
        return methods
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _defineValues(self):
        """ This function get some parameters of the micrographs"""
        # Get pointer to input micrographs 
        self.inputMics = self.inputMicrographs.get() 
        acquisition = self.inputMics.getAcquisition()
        
        self._params = {'voltage': acquisition.getVoltage(),
                        'sphericalAberration': acquisition.getSphericalAberration(),
                        'magnification': acquisition.getMagnification(),
                        'ampContrast': acquisition.getAmplitudeContrast(),
                        'samplingRate': self.inputMics.getSamplingRate(),
                        'scannedPixelSize': self.inputMics.getScannedPixelSize(),
                        'windowSize': self.windowSize.get(),
                        'lowRes': self.lowRes.get(),
                        'highRes': self.highRes.get(),
                        # Convert from microns to Amstrongs
                        'minDefocus': self.minDefocus.get() * 1e+4, 
                        'maxDefocus': self.maxDefocus.get() * 1e+4
                       }
    
    def _defineRecalValues(self, ctfModel):
        """ This function get the acquisition info of the micrographs"""
        mic = ctfModel.getMicrograph()
        
        acquisition = mic.getAcquisition()
        scannedPixelSize = mic.getSamplingRate() * acquisition.getMagnification() / 10000
        self._params = {'voltage': acquisition.getVoltage(),
                        'sphericalAberration': acquisition.getSphericalAberration(),
                        'magnification': acquisition.getMagnification(),
                        'ampContrast': acquisition.getAmplitudeContrast(),
                        'scannedPixelSize': scannedPixelSize,
                        'samplingRate': mic.getSamplingRate()
                       }
    
    def _getPrevMicDir(self, ctfModel):
        return dirname(ctfModel.getPsdFile())
    
    def _ctfCounter(self, values):
        """ This function return the number of CTFs that was recalculated.
        """
        numberOfCTF = len(values)/2
        msg = "CTF Re-estimation of %d micrographs" % numberOfCTF
        self.summaryVar.set(msg)
    
    def _getInputCtf(self):
        if self.continueRecal:
            sqliteFile = self._getPath()
#             return self.outputCTF.get()
        else:
            return self.inputCtf.get()
        
    def _getMicrographDir(self, mic):
        """ Return an unique dir name for results of the micrograph. """
        return self._getExtraPath(removeBaseExt(mic.getFileName()))        
    
    def _iterMicrographs(self):
        """ Iterate over micrographs and yield
        micrograph name and a directory to process.
        """
        for mic in self.inputMics:
            micFn = mic.getFileName()
            micDir = self._getMicrographDir(mic)
            yield (micFn, micDir, mic)  
    
    def _prepareCommand(self):
        """ This function should be implemented to prepare the
        arguments template if doesn't change for each micrograph
        After this method self._program and self._args should be set. 
        """
        pass
    
    def _defocusMaxMin(self, defocusList):
        """ This function return the minimum and maximum of the defocus
        of a SetOfMicrographs.
        """
        minimum = float(min(defocusList))/10000
        maximum = float(max(defocusList))/10000
        msg = "The range of micrograph's experimental defocus was %(minimum)0.3f - %(maximum)0.3f microns. " % locals()

        self.methodsVar.set(msg)


class ProtPreprocessMicrographs(ProtMicrographs):
    pass

