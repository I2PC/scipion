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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from os.path import exists, dirname, join

from pyworkflow.object import Boolean
from pyworkflow.protocol.constants import (STEPS_PARALLEL, LEVEL_ADVANCED,
                                           STATUS_NEW)
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,
                                        BooleanParam, FileParam)
from pyworkflow.utils.path import copyTree, removeBaseExt, makePath
from pyworkflow.utils.properties import Message
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.data import SetOfMicrographs, SetOfCTF



class ProtMicrographs(EMProtocol):
    pass


class ProtCTFMicrographs(ProtMicrographs):
    """ Base class for all protocols that estimates the CTF"""

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL
        self.isFirstTime = Boolean(False)

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_CTF_ESTI)
        form.addParam('recalculate', BooleanParam, default=False,
                      condition='recalculate',
                      label="Do recalculate ctf?")

        form.addParam('continueRun', PointerParam, allowsNull=True,
                      condition='recalculate', label="Input previous run",
                      pointerClass=self.getClassName())
        form.addHidden('sqliteFile', FileParam, condition='recalculate',
                       allowsNull=True)

        form.addParam('inputMicrographs', PointerParam, important=True,
                      condition='not recalculate', label=Message.LABEL_INPUT_MIC,
                      pointerClass='SetOfMicrographs')
        form.addParam('ctfDownFactor', FloatParam, default=1.,
                      label='CTF Downsampling factor',
                      condition='not recalculate',
                      help='Set to 1 for no downsampling. Non-integer downsample '
                           'factors are possible. This downsampling is only used '
                           'for estimating the CTF and it does not affect any '
                           'further calculation. Ideally the estimation of the '
                           'CTF is optimal when the Thon rings are not too '
                           'concentrated at the origin (too small to be seen) '
                           'and not occupying the whole power spectrum (since '
                           'this downsampling might entail aliasing).')

        self._defineProcessParams(form)

        line = form.addLine('Resolution', condition='not recalculate',
                            help='Give a value in digital frequency '
                                 '(i.e. between 0.0 and 0.5). These cut-offs '
                                 'prevent the typical peak at the center of the'
                                 ' PSD and high-resolution terms where only '
                                 'noise exists, to interfere with CTF '
                                 'estimation. The default lowest value is 0.05 '
                                 'but for micrographs with a very fine sampling '
                                 'this may be lowered towards 0. The default '
                                 'highest value is 0.35, but it should be '
                                 'increased for micrographs with signals '
                                 'extending beyond this value. However, if '
                                 'your micrographs extend further than 0.35, '
                                 'you should consider sampling them at a finer '
                                 'rate.')
        line.addParam('lowRes', FloatParam, default=0.05, label='Lowest' )
        line.addParam('highRes', FloatParam, default=0.35, label='Highest')
        line = form.addLine('Defocus search range (microns)',
                            condition='not recalculate',
                            expertLevel=LEVEL_ADVANCED,
                            help='Select _minimum_ and _maximum_ values for '
                                 'defocus search range (in microns). Underfocus'
                                 ' is represented by a positive number.')
        line.addParam('minDefocus', FloatParam, default=0.25,
                      label='Min')
        line.addParam('maxDefocus', FloatParam, default=4.,
                      label='Max')

        form.addParam('windowSize', IntParam, default=256,
                      expertLevel=LEVEL_ADVANCED,
                      label='Window size', condition='not recalculate',
                      help='The PSD is estimated from small patches of this '
                           'size. Bigger patches allow identifying more '
                           'details. However, since there are fewer windows, '
                           'estimations are noisier.')

        form.addParallelSection(threads=2, mpi=1)

    def _defineProcessParams(self, form):
        """ This method should be implemented by subclasses
        to add other parameter relatives to the specific operation."""
        pass

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        """ Insert the steps to perform CTF estimation, or re-estimation,
        on a set of micrographs.
        """
        # Store all steps ids, final step createOutput depends on all of them
        deps = []
        fDeps = []
        self.insertedDict = {}

        if not self.recalculate:
            deps = self._insertEstimationSteps(self.insertedDict,
                                               self.inputMicrographs.get())
            # Insert step to create output objects
            fDeps = self._insertFinalSteps(deps)
            # For the streaming mode, the steps function have a 'wait' flag
            # that can be turned on/off. For example, here we insert the
            # createOutputStep but it wait=True, which means that can not be
            # executed until it is set to False
            # (when the input micrographs stream is closed)
            waitCondition = self._getFirstJoinStepName() == 'createOutputStep'
        else:
            if self.isFirstTime:
                # Insert previous estimation or re-estimation an so on...
                self._insertPreviousSteps()
                self.isFirstTime.set(False)
            fDeps = self._insertRecalculateSteps()
            # For now the streaming is not allowed for recalculate CTF
            waitCondition = False

        self._insertFunctionStep('createOutputStep', prerequisites=fDeps,
                                 wait=waitCondition)

    def _insertFinalSteps(self, deps):
        """ This should be implemented in subclasses"""
        return deps

    def _getFirstJoinStepName(self):
        # This function will be used for streamming, to check which is
        # the first function that need to wait for all micrographs
        # to have completed, this can be overriden in subclasses
        # (e.g., in Xmipp 'sortPSDStep')
        return 'createOutputStep'

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None

    def _checkNewMicrographs(self, micSet, outputStep):
        """ Check if there are new micrographs to be processed
        and add the necessary steps.
        """
        newMics = []
        for micFn, _, mic in self._iterMicrographs(micSet):
            if mic.getMicName() not in self.insertedDict:
                newMics.append(micFn)

        if newMics:
            fDeps = self._insertEstimationSteps(self.insertedDict, micSet)
            self._storeSteps()
            self._numberOfSteps.set(len(self._steps))
            self._store(self._numberOfSteps)
            if outputStep:
                outputStep.addPrerequisites(*fDeps)

        return newMics

    def _checkNewCTFs(self, micSet):
        """ Check for already computed CTF and update the output set. """
        newCTFs = []
        ctfDict = {}
        ctfSet = SetOfCTF(filename=self._getPath('ctfs.sqlite'))
        for ctf in ctfSet:
            ctfDict[ctf.getObjId()] = True

        if ctfDict: # it means there are previous ctfs computed
            if ctfSet.getSize():
                ctfSet.enableAppend()
        else:
            ctfSet.setStreamState(ctfSet.STREAM_OPEN)

        for micFn, micDir, mic in self._iterMicrographs(micSet):
            if (exists(self._getMicrographDone(micDir))
                and not mic.getObjId() in ctfDict):
                ctf = self._createCtfModel(mic)
                ctfSet.append(ctf)
                newCTFs.append(mic.getObjId())

        return ctfSet, newCTFs

    def _stepsCheck(self):
        # For now the streaming is not allowed for recalculate CTF
        if self.recalculate:
            return

        # Check if there are new micrographs to process
        micFn = self.inputMicrographs.get().getFileName()
        micSet = SetOfMicrographs(filename=micFn)
        micSet.loadAllProperties()
        streamClosed = micSet.isStreamClosed()
        endCTFs = False

        outputStep = self._getFirstJoinStep()
        self._checkNewMicrographs(micSet, outputStep)
        ctfSet, newCTFs = self._checkNewCTFs(micSet)

        endCTFs = streamClosed and micSet.getSize() == ctfSet.getSize()

        if newCTFs:
            if endCTFs:
                ctfSet.setStreamState(ctfSet.STREAM_CLOSED)
            self._updateOutput(ctfSet)

        if outputStep and outputStep.isWaiting() and endCTFs:
            outputStep.setStatus(STATUS_NEW)

        micSet.close()

    def _insertEstimationSteps(self, insertedDict, inputMics):
        estimDeps = []
        self._defineValues()
        self._prepareCommand()
        # For each micrograph insert the steps to process it
        for micFn, micDir, mic in self._iterMicrographs(inputMics):
            if mic.getMicName() not in insertedDict:
                # CTF estimation
                # Make estimation steps independent between them
                stepId = self._insertFunctionStep('_estimateCTF', micFn, micDir,
                                                  mic.getMicName(),
                                                  prerequisites=[])
                estimDeps.append(stepId)
                insertedDict[mic.getMicName()] = stepId
        return estimDeps

    def _insertRecalculateSteps(self):
        recalDeps = []
        # For each psd insert the steps to process it
        self.recalculateSet = SetOfCTF(filename=self.sqliteFile.get(),
                                       objDoStore=False)
        for ctf in self.recalculateSet:
            line = ctf.getObjComment()
            if ctf.isEnabled() and line:
                # CTF Re-estimation
                copyId = self._insertFunctionStep('copyMicDirectoryStep',
                                                  ctf.getObjId())
                # Make estimation steps independent between them
                stepId = self._insertFunctionStep('_restimateCTF',
                                                  ctf.getObjId(),
                                                  prerequisites=[copyId])
                recalDeps.append(stepId)
        return recalDeps

    #--------------------------- STEPS functions -------------------------------
    def _estimateCTF(self, micFn, micDir, micName):
        """ Do the CTF estimation with the specific program
        and the parameters required.
        Params:
         micFn: micrograph filename
         micDir: micrograph directory
        """
        raise Exception(Message.ERROR_NO_EST_CTF)

    def _restimateCTF(self, micId):
        """ Do the CTF estimation with the specific program
        and the parameters required.
        Params:
         micFn: micrograph filename
         micDir: micrograph directory
        """
        raise Exception(Message.ERROR_NO_EST_CTF)

    def copyMicDirectoryStep(self, micId):
        """ Copy micrograph's directory tree for recalculation"""
        ctfModel = self.recalculateSet[micId]
        mic = ctfModel.getMicrograph()

        prevDir = self._getPrevMicDir(ctfModel)
        micDir = self._getMicrographDir(mic)
        if not prevDir == micDir:
            # Create micrograph dir under extra directory
            makePath(micDir)
            if not exists(micDir):
                raise Exception("No created dir: %s " % micDir)
            copyTree(prevDir, micDir)

    def _createCtfModel(self, mic):
        """ This should be implemented in subclasses
        in order to create a CTF model from program results.
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
            prot = self.continueRun.get() or self
            micSet = prot.outputCTF.getMicrographs()
            # We suppose this is reading the ctf selection
            # (with enabled/disabled) to only consider the enabled ones
            # in the final SetOfCTF
            #TODO: maybe we can remove the need of the extra text file
            # with the recalculate parameters
            newCount = 0
            for ctfModel in self.recalculateSet:
                if ctfModel.isEnabled() and ctfModel.getObjComment():
                    mic = ctfModel.getMicrograph()
                    # Update the CTF models that where recalculated and append
                    # later to the set, we don't want to copy the id here since
                    # it is already correct
                    newCtf = self._createCtfModel(mic, updateSampling=False)
                    ctfModel.copy(newCtf, copyId=False)
                    ctfModel.setEnabled(True)
                    newCount += 1
                ctfSet.append(ctfModel)
            ctfSet.setMicrographs(micSet)
            self._defineOutputs(outputCTF=ctfSet)
            self._defineCtfRelation(micSet, ctfSet)
            self._computeDefocusRange(ctfSet)
            self.summaryVar.set("CTF Re-estimation of %d micrographs"
                                % newCount)
        else:
            self._createOutputStep()

    #--------------------------- INFO functions --------------------------------
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
                summary.append("CTF estimation of %d micrographs."
                               % self.inputMicrographs.get().getSize())

        return summary

    def _methods(self):
        methods = []

        if hasattr(self, 'outputCTF') and self.isFinished():
            methods.append(self.methodsVar.get())
        else:
            methods.append(Message.TEXT_NO_CTF_READY)

        return methods

    #--------------------------- UTILS functions -------------------------------
    def _defineValues(self):
        """ This function get some parameters of the micrographs"""
        # Get pointer to input micrographs 
        self.inputMics = self.inputMicrographs.get()
        acq = self.inputMics.getAcquisition()

        self._params = {'voltage': acq.getVoltage(),
                        'sphericalAberration': acq.getSphericalAberration(),
                        'magnification': acq.getMagnification(),
                        'ampContrast': acq.getAmplitudeContrast(),
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

        acq = mic.getAcquisition()
        mag = acq.getMagnification()
        scannedPixelSize = mic.getSamplingRate() * mag / 10000
        self._params = {'voltage': acq.getVoltage(),
                        'sphericalAberration': acq.getSphericalAberration(),
                        'magnification': mag,
                        'ampContrast': acq.getAmplitudeContrast(),
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

    def _getMicrographDone(self, micDir):
        """ Return the file that is used as a flag of termination. """
        return join(micDir, 'done.txt')

    def _writeMicrographDone(self, micDir):
        open(self._getMicrographDone(micDir), 'w').close()

    def _iterMicrographs(self, inputMics=None):
        """ Iterate over micrographs and yield
        micrograph name and a directory to process.
        """
        if inputMics is None:
            inputMics = self.inputMics

        for mic in inputMics:
            micFn = mic.getFileName()
            micDir = self._getMicrographDir(mic)
            yield (micFn, micDir, mic)

    def _prepareCommand(self):
        """ This function should be implemented to prepare the
        arguments template if doesn't change for each micrograph
        After this method self._program and self._args should be set. 
        """
        pass

    def _computeDefocusRange(self, ctfSet):
        """ Compute the minimum and maximu defoucs in a set of CTFs.
        The protocol methodsVar will be updated with new values.

        Params:
            ctfSet: the set of CTFs to compute min and max
        """
        defocusList = []

        for ctf in ctfSet:
            defocusList.append(ctf.getDefocusU())
            defocusList.append(ctf.getDefocusV())

        minD = min(defocusList) / 10000.
        maxD = max(defocusList) / 10000.

        self.methodsVar.set("The range of micrograph's experimental defocus was"
                            " %0.3f - %0.3f microns. " % (minD, maxD))

        self._store(self.methodsVar)

    def _defocusMaxMin(self, defocusList):
        """ This function return the minimum and maximum of the defocus
        of a SetOfMicrographs.
        """
        raise Exception("DEPRECATED")


class ProtPreprocessMicrographs(ProtMicrographs):
    pass

