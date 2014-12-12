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

from pyworkflow.object import String
from pyworkflow.protocol.constants import STEPS_PARALLEL, LEVEL_ADVANCED, LEVEL_EXPERT
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
        self.methodsInfo = String()
        self.stepsExecutionMode = STEPS_PARALLEL
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_CTF_ESTI)
        
        form.addParam('inputMicrographs', PointerParam, important=True,
                      label=Message.LABEL_INPUT_MIC, pointerClass='SetOfMicrographs')
        form.addParam('ctfDownFactor', FloatParam, default=1.,
                      label='CTF Downsampling factor',
                      help='Set to 1 for no downsampling. Non-integer downsample factors are possible. '
                      'This downsampling is only used for estimating the CTF and it does not affect '
                      'any further calculation. Ideally the estimation of the CTF is optimal when '
                      'the Thon rings are not too concentrated at the origin (too small to be seen) '
                      'and not occupying the whole power spectrum (since this downsampling might '
                      'entail aliasing).')
        
        self._defineProcessParams(form)

        line = form.addLine('Resolution', 
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
                            help='Select _minimum_ and _maximum_ values for defocus search range (in microns).'
                                 'Underfocus is represented by a positive number.')
        line.addParam('minDefocus', FloatParam, default=0.5, 
                      label='Min')
        line.addParam('maxDefocus', FloatParam, default=4.,
                      label='Max')
        
        form.addParam('windowSize', IntParam, default=256, expertLevel=LEVEL_ADVANCED,
                      label='Window size', 
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
        """ Insert the steps to perform ctf estimation on a set of micrographs.
        """
        self._defineValues()
        self._prepareCommand()
        deps = [] # Store all steps ids, final step createOutput depends on all of them
        # For each micrograph insert the steps to process it
        for micFn, micDir, _ in self._iterMicrographs():
            # CTF estimation
            # Make estimation steps independent between them
            stepId = self._insertFunctionStep('_estimateCTF', micFn, micDir,
                                                  prerequisites=[]) # Make estimation steps independent between them
            deps.append(stepId)
        self._insertFinalSteps(deps)
    
    
    def _insertFinalSteps(self, deps):
        # Insert step to create output objects       
        self._insertFunctionStep('createOutputStep', prerequisites=deps)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _estimateCTF(self, micFn, micDir):
        """ Do the CTF estimation with the specific program
        and the parameters required.
        Params:
         micFn: micrograph filename
         micDir: micrograph directory
        """
        raise Exception(Message.ERROR_NO_EST_CTF)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputCTF'):
            summary.append(Message.TEXT_NO_CTF_READY)
        else:
            summary.append("CTF estimation of %d micrographs." % self.inputMicrographs.get().getSize())
        return summary
    
    def _methods(self):
        methods = []
        
        if not hasattr(self, 'outputCTF'):
            methods.append(Message.TEXT_NO_CTF_READY)
        else:
            methods.append(self.methodsInfo.get())
            
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
        self.methodsInfo.set(msg)


class ProtRecalculateCTF(ProtMicrographs):
    """ Base class for all protocols that re-calculate the CTF"""
    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.methodsInfo = String()
        self.summaryInfo = String()
        self.stepsExecutionMode = STEPS_PARALLEL
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_CTF_ESTI)
        
        form.addParam('inputCtf', PointerParam, important=True,
              label="input the SetOfCTF to recalculate", pointerClass='SetOfCTF')
        form.addHidden('sqliteFile', FileParam)
        
        form.addParallelSection(threads=1, mpi=1)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ Insert the steps to perform ctf re-estimation on a set of CTFs.
        """
        #self._insertFunctionStep('createSubsetOfCTF')
        deps = [] # Store all steps ids, final step createOutput depends on all of them
        # For each psd insert the steps to process it
        self.recalculateSet = SetOfCTF(filename=self.sqliteFile.get(), objDoStore=False)
        for ctf in self.recalculateSet:
            line = ctf.getObjComment()
            if ctf.isEnabled() and line:
                # CTF Re-estimation
                copyId = self._insertFunctionStep('copyMicDirectory', ctf.getObjId())
                stepId = self._insertFunctionStep('_estimateCTF', ctf.getObjId(), prerequisites=[copyId]) # Make estimation steps independent between them
                deps.append(stepId)
        # Insert step to create output objects
        self._insertFinalSteps(deps)


    def _insertFinalSteps(self, deps):
        # Insert step to create output objects       
        self._insertFunctionStep('createOutputStep', prerequisites=deps)

    #--------------------------- STEPS functions ---------------------------------------------------

    def copyMicDirectory(self, id):
        """ Copy micrograph's directory tree for recalculation"""
        ctfModel = self.recalculateSet[id]
        mic = ctfModel.getMicrograph()
        
        prevDir = self._getPrevMicDir(ctfModel)
        micDir = self._getMicrographDir(mic)
        # Create micrograph dir under extra directory
        makePath(micDir)
        if not exists(micDir):
            raise Exception("No created dir: %s " % micDir)
        copyTree(prevDir, micDir)
    
    def _estimateCTF(self, id):
        """ Do the CTF estimation with the specific program
        and the parameters required.
        Params:
         micFn: micrograph filename
         micDir: micrograph directory
        """
        raise Exception(Message.ERROR_NO_EST_CTF)
    
    def _createNewCtfModel(self, mic):
        """ This should be implemented in subclasses
        in order to create a CTF model 
        """
        pass
    
    def createOutputStep(self):
        """ This function is shared by Xmipp and CTFfind recalculate
        protocols, it will iterated for each CTF model, see
        if was recalculated and update with new defocus values.
        The function that should be implemented in each subclass
        is _createOutputStep.
        """
        ctfSet = self._createSetOfCTF("_recalculated")
        defocusList = []
        # README: We suppose this is reading the ctf selection (with enabled/disabled)
        # to only consider the enabled ones in the final SetOfCTF

        
        #TODO: maybe we can remove the need of the extra text file
        # with the recalculate parameters
        for ctfModel, ctfModel2 in izip(self.inputCtf.get(), self.recalculateSet):
            if ctfModel2.isEnabled() and ctfModel2.getObjComment():
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
            
        self._defineOutputs(outputCTF=ctfSet)
        self._defineSourceRelation(self.inputCtf.get(), ctfSet)

        self._defocusMaxMin(defocusList)
        self._ctfCounter(defocusList)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _summary(self):
        summary = []
        if not (hasattr(self, 'outputCTF') ):
            summary.append(Message.TEXT_NO_CTF_READY)
        else:
            summary.append(self.summaryInfo.get())
        return summary
    
    def _methods(self):
        methods = []
        
        if not (hasattr(self, 'outputCTF') ):
            methods.append(Message.TEXT_NO_CTF_READY)
        else:
            methods.append(self.methodsInfo.get())
            
        return methods
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _defineValues(self, ctfModel):
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
    
    def _getMicrographDir(self, mic):
        """ Return an unique dir name for results of the micrograph. """
        return self._getExtraPath(removeBaseExt(mic.getFileName()))        
    
    def _prepareCommand(self):
        """ This function should be implemented to prepare the
        arguments template if doesn't change for each micrograph
        After this method self._program and self._args should be set. 
        """
        pass
    
    def _getPrevMicDir(self, ctfModel):
        return dirname(ctfModel.getPsdFile())
    


    def _defocusMaxMin(self, values):
        """ This function return the minimum and maximum of the defocus
        of a SetOfMicrographs.
        """
        minimum = float(min(values))/10000
        maximum = float(max(values))/10000
        msg = "The range of micrograph's experimental defocus was %(minimum)0.3f - %(maximum)0.3f microns. " % locals()
        self.methodsInfo.set(msg)
    
    def _ctfCounter(self, values):
        """ This function return the number of CTFs that was recalculated.
        """
        numberOfCTF = len(values)/2
        msg = "CTF Re-estimation of %(numberOfCTF)d micrographs" % locals()
        self.summaryInfo.set(msg)


class ProtPreprocessMicrographs(ProtMicrographs):
    pass

