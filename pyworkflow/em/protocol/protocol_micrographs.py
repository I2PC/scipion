# **************************************************************************
# *
# * Authors:     Josue Gomez BLanco (jgomez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *              Airen Zaldivar Peraza (azaldivar@cnb.csic.es)
# *
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

from os.path import exists, dirname, join, getmtime
from datetime import datetime
from collections import OrderedDict


from pyworkflow.object import Set, Boolean, Pointer
from pyworkflow.protocol.constants import (STEPS_PARALLEL, LEVEL_ADVANCED,
                                           STATUS_NEW)
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,
                                        BooleanParam, FileParam)
from pyworkflow.utils.path import copyTree, removeBaseExt, makePath, makeFilePath
from pyworkflow.utils.properties import Message
from pyworkflow.utils.utils import prettyTime
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
                      condition='not recalculate',
                      label=Message.LABEL_INPUT_MIC,
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
        if not self.recalculate:
            self.initialIds = self._insertInitialSteps()
            self.micDict = OrderedDict()
            micDict, _ = self._loadInputList()
            ctfIds = self._insertNewMicsSteps(micDict.values())
            self._insertFinalSteps(ctfIds)
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
            ctfIds = self._insertRecalculateSteps()
            # For now the streaming is not allowed for recalculate CTF
            waitCondition = False

        self._insertFunctionStep('createOutputStep', prerequisites=ctfIds,
                                 wait=waitCondition)

    def _insertInitialSteps(self):
        """ Override this function to insert some steps before the
        estimate ctfs steps.
        Should return a list of ids of the initial steps. """
        return []

    def _insertNewMicsSteps(self, inputMics):
        """ Insert steps to process new mics (from streaming)
        Params:
            inputMics: input mics set to be check
        """
        deps = []
        # For each mic insert the step to process it
        for mic in inputMics:
            micKey = mic.getMicName()
            if micKey not in self.micDict:
                args = [mic.getFileName(), self._getMicrographDir(mic), micKey]
                stepId = self._insertEstimationSteps(self.initialIds, *args)
                deps.append(stepId)
                self.micDict[micKey] = mic
        return deps

    def _insertEstimationSteps(self, prerequisites, *args):
        """ Basic method to insert a estimateCTF step for a given micrograph."""
        self._defineValues()
        self._prepareCommand()
        micStepId = self._insertFunctionStep('_estimateCTF', *args,
                                             prerequisites=prerequisites)
        return micStepId

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
        self.inputMics = self.getInputMicrographs()
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

        self.methodsVar.set("Estimated  defocus range defocus was"
                            " %0.3f - %0.3f microns. " % (minD, maxD))

        self._store(self.methodsVar)

    def _defocusMaxMin(self, defocusList):
        """ This function return the minimum and maximum of the defocus
        of a SetOfMicrographs.
        """
        raise Exception("DEPRECATED")

    def getInputMicrographsPointer(self):
        return self.inputMicrographs

    def getInputMicrographs(self):
        return self.getInputMicrographsPointer().get()

    # ------ Methods for Streaming picking --------------
    def _stepsCheck(self):
        # To allow streaming ctf estimation we need to detect:
        #   1) new micrographs ready to be picked
        #   2) new output ctfs that have been produced and add then
        #      to the output set.
    
        # For now the streaming is not allowed for recalculate CTF
        if self.recalculate:
            return
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        # Check if there are new micrographs to process from the input set
        localFile = self.getInputMicrographs().getFileName()
        now = datetime.now()
        self.lastCheck = getattr(self, 'lastCheck', now)
        mTime = datetime.fromtimestamp(getmtime(localFile))
        self.debug('Last check: %s, modification: %s'
                  % (prettyTime(self.lastCheck),
                     prettyTime(mTime)))
        # If the input micrographs.sqlite have not changed since our last check,
        # it does not make sense to check for new input data
        if self.lastCheck > mTime and hasattr(self, 'listOfMics'):
            return None

        self.lastCheck = now
        # Open input micrographs.sqlite and close it as soon as possible
        micDict, self.streamClosed = self._loadInputList()
        newMics = micDict.values()
        outputStep = self._getFirstJoinStep()

        if newMics:
            fDeps = self._insertNewMicsSteps(newMics)
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()

    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return
        # Load previously done items (from text file)
        doneList = self._readDoneList()
        # Check for newly done items
        listOfMics = self.micDict.values()
        nMics = len(listOfMics)
        newDone = [m for m in listOfMics
                   if m.getObjId() not in doneList and self._isMicDone(m)]

        # Update the file with the newly done mics
        # or exit from the function if no new done mics
        self.debug('_checkNewOutput: ')
        self.debug('   listOfMics: %s, doneList: %s, newDone: %s'
                   % (nMics, len(doneList), len(newDone)))

        allDone = len(doneList) + len(newDone)
        # We have finished when there is not more input mics (stream closed)
        # and the number of processed mics is equal to the number of inputs
        self.finished = self.streamClosed and allDone == nMics
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN
        self.debug('   streamMode: %s newDone: %s' % (streamMode,
                                                      not(newDone == [])))

        if newDone:
            newDoneUpdated = self._updateOutputCTFSet(newDone, streamMode)
            self._writeDoneList(newDoneUpdated)
        elif not self.finished:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        self.debug('   finished: %s ' % self.finished)
        self.debug('        self.streamClosed (%s) AND' % self.streamClosed)
        self.debug('        allDone (%s) == len(self.listOfMics (%s)'
                   % (allDone, nMics))

        if self.finished:  # Unlock createOutputStep if finished all jobs
            self._updateStreamState(streamMode)
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)

    def _loadInputList(self):
        """ Load the input set of micrographs that are ready to be picked. """
        return self._loadSet(self.getInputMicrographs(), SetOfMicrographs,
                        lambda mic: mic.getMicName())
        
    def _loadSet(self, inputSet, SetClass, getKeyFunc):
        """ Load a given input set if their items are not already present
        in the self.micDict.
        This can be used to load new micrographs for picking as well as
        new CTF (if used) in streaming.
        """
        setFn = inputSet.getFileName()
        self.debug("Loading input db: %s" % setFn)
        updatedSet = SetClass(filename=setFn)
        updatedSet.loadAllProperties()
        newItemDict = OrderedDict()
        for item in updatedSet:
            micKey = getKeyFunc(item)
            if micKey not in self.micDict:
                newItemDict[micKey] = item.clone()
        streamClosed = updatedSet.isStreamClosed()
        updatedSet.close()
        self.debug("Closed db.")

        return newItemDict, streamClosed

    def _updateOutputCTFSet(self, micList, streamMode):
        micDoneList = [mic for mic in micList]
        # Do no proceed if there is not micrograph ready
        if not micDoneList:
            return []

        outputName = 'outputCTF'
        outputCtf = getattr(self, outputName, None)

        # If there is not outputCTF yet, it means that is the first
        # time we are updating output CTFs, so we need to first create
        # the output set
        firstTime = outputCtf is None

        if firstTime:
            outputCtf = self._createSetOfCTF()
        else:
            outputCtf.enableAppend()

        outputCtf.setMicrographs(self.getInputMicrographs())

        for micFn, micDir, mic in self._iterMicrographs(micList):
            ctf = self._createCtfModel(mic)
            outputCtf.append(ctf)

        self.debug(" _updateOutputCTFSet Stream Mode: %s " % streamMode)
        self._updateOutputSet(outputName, outputCtf, streamMode)

        if firstTime:  # define relation just once
            # Using a pointer to define the relations is more robust to
            # scheduling and id changes between the protocol run.db and
            # the main project database.
            self._defineCtfRelation(self.getInputMicrographsPointer(),
                                    outputCtf)

        return micDoneList

    def _updateStreamState(self, streamMode):
        outputName = 'outputCTF'
        outputCtf = getattr(self, outputName, None)

        # If there are not outputCoordinates yet, it means that is the first
        # time we are updating output coordinates, so we need to first create
        # the output set
        firstTime = outputCtf is None

        if firstTime:
            micSetPtr = self.getInputMicrographsPointer()
            outputCtf = self._createSetOfCoordinates(micSetPtr)
        else:
            outputCtf.enableAppend()

        self.debug(" _updateStreamState Stream Mode: %s " % streamMode)
        self._updateOutputSet(outputName, outputCtf, streamMode)

    def _readDoneList(self):
        """ Read from a text file the id's of the items that have been done. """
        doneFile = self._getAllDone()
        doneList = []
        # Check what items have been previously done
        if exists(doneFile):
            with open(doneFile) as f:
                doneList += [int(line.strip()) for line in f]
        return doneList

    def _writeDoneList(self, micList):
        """ Write to a text file the items that have been done. """
        doneFile = self._getAllDone()

        if not exists(doneFile):
            makeFilePath(doneFile)

        with open(doneFile, 'a') as f:
            for mic in micList:
                f.write('%d\n' % mic.getObjId())

    def _isMicDone(self, mic):
        """ A mic is done if the marker file exists. """
        micDir = self._getMicrographDir(mic)
        return exists(self._getMicrographDone(micDir))

    def _getAllDone(self):
        return self._getExtraPath('DONE', 'all.TXT')


class ProtPreprocessMicrographs(ProtMicrographs):
    pass

