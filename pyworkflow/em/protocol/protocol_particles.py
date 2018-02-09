# **************************************************************************
# *
# * Authors:     Airen Zaldivar Peraza (azaldivar@cnb.csic.es) [1]
# *              J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] SciLifeLab, Stockholm University
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

import os
from datetime import datetime
from collections import OrderedDict

from pyworkflow.object import Set, String, Pointer
import pyworkflow.protocol.params as params
from pyworkflow.protocol import STATUS_NEW
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.constants import RELATION_CTF
from pyworkflow.em.data import (EMObject, SetOfCoordinates, Micrograph,
                                SetOfMicrographs, SetOfCTF)

import pyworkflow.utils as pwutils
from pyworkflow.utils.properties import Message



class ProtParticles(EMProtocol):
    pass


class ProtProcessParticles(ProtParticles):
    """ This class will serve as a base for all protocol
    that performs some operation on Particles (i.e. filters, mask, resize, etc)
    It is mainly defined by an inputParticles and outputParticles.
    """
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        
        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      label=Message.LABEL_INPUT_PART, important=True)
        # Hook that should be implemented in subclasses
        self._defineProcessParams(form)
        
        __threads, __mpi = self._getDefaultParallel()
        
        form.addParallelSection(threads=__threads, mpi=__mpi)
        
    def _defineProcessParams(self, form):
        """ This method should be implemented by subclasses
        to add other parameter relatives to the specific operation."""
        pass
    
    def _getDefaultParallel(self):
        """ Return the default value for thread and MPI
        for the parallel section definition.
        """
        return (0, 0)


class ProtFilterParticles(ProtProcessParticles):
    """ Base class for filters on particles of type ProtPreprocessParticles.
    """
    pass


class ProtOperateParticles(ProtProcessParticles):
    """ Base class for operations on particles of type ProtPreprocessParticles.
    """
    def __init__(self, **args):
        ProtProcessParticles.__init__(self, **args)


class ProtMaskParticles(ProtProcessParticles):
    """ This is the base for the branch of mask, 
    between the ProtPreprocessParticles """
    pass


class ProtParticlePicking(ProtParticles):
    OUTPUT_PREFIX = 'outputCoordinates'

    def _defineParams(self, form):

        form.addSection(label='Input')
        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      label=Message.LABEL_INPUT_MIC, important=True,
                      help='Select the SetOfMicrographs to be used during '
                           'picking.')

    #--------------------------- INFO functions --------------------------------
    def getSummary(self, coordSet):
        summary = []
        summary.append("Number of particles picked: %s" % coordSet.getSize())
        summary.append("Particle size: %s" % coordSet.getBoxSize())
        return "\n".join(summary)

    def getMethods(self, output):
        msg = 'User picked %d particles ' % output.getSize()
        msg += 'with a particle size of %s.' % output.getBoxSize()
        return msg

    def _methods(self):
        methodsMsgs = []
        if self.getInputMicrographs() is None:
            return ['Input micrographs not available yet.']

        methodsMsgs.append("Input micrographs %s of size %d."
                           % (self.getObjectTag(self.getInputMicrographs()),
                              self.getInputMicrographs().getSize()))

        if self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes(EMObject):
                msg = self.getMethods(output)
                methodsMsgs.append("%s: %s"%(self.getObjectTag(output), msg))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs

    def _summary(self):
        summary = []
        if self.getInputMicrographs() is not None:
            summary.append("Number of input micrographs: %d"
                           % self.getInputMicrographs().getSize())

        if self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes(EMObject):
                summary.append("*%s:* \n %s " % (key, output.getObjComment()))
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary

    def getInputMicrographsPointer(self):
        return self.inputMicrographs

    def getInputMicrographs(self):
        return self.getInputMicrographsPointer().get()

    def _getCoords(self, CoordClass):
        result = None
        for _, attr in self.iterOutputAttributes(CoordClass):
            result = attr # Get the last output that is SetOfCoordinates or so
        return result

    def getCoords(self):
        return self._getCoords(SetOfCoordinates)

    def getCoordsTiltPair(self):
        from pyworkflow.em.data_tiltpairs import CoordinatesTiltPair
        return self._getCoords(CoordinatesTiltPair)

    def _createOutput(self, outputDir):
        micSet = self.getInputMicrographs()
        suffix = self.__getOutputSuffix()
        outputName = self.OUTPUT_PREFIX + suffix
        coordSet = self._createSetOfCoordinates(micSet, suffix)
        self.readSetOfCoordinates(outputDir, coordSet)
        coordSet.setObjComment(self.getSummary(coordSet))
        outputs = {outputName: coordSet}
        self._defineOutputs(**outputs)
        self._defineSourceRelation(self.getInputMicrographsPointer(), coordSet)

    def createOutputStep(self):
        self._createOutput(self._getExtraPath())

    def readSetOfCoordinates(self, workingDir, coordSet):
        pass

    def getCoordsDir(self):
        pass

    def __getOutputSuffix(self):
        """ Get the name to be used for a new output.
        For example: outputCoordiantes7.
        It should take into account previous outputs
        and number with a higher value.
        """
        maxCounter = -1
        for attrName, _ in self.iterOutputAttributes(SetOfCoordinates):
            suffix = attrName.replace(self.OUTPUT_PREFIX, '')
            try:
                counter = int(suffix)
            except:
                counter = 1 # when there is not number assume 1
            maxCounter = max(counter, maxCounter)

        return str(maxCounter+1) if maxCounter > 0 else '' # empty if not output

    def registerCoords(self, coordsDir):
        """ This method is usually inherited by all Pickers
        and it is used from the Java picking GUI to register
        a new SetOfCoordinates when the user click on +Particles button. 
        """
        suffix = self.__getOutputSuffix()
        outputName = self.OUTPUT_PREFIX + suffix

        from pyworkflow.em.packages.xmipp3 import readSetOfCoordinates
        inputset = self.getInputMicrographs()
        # micrographs are the input set if protocol is not finished
        outputset = self._createSetOfCoordinates(inputset, suffix=suffix)
        readSetOfCoordinates(coordsDir, outputset.getMicrographs(), outputset)
        summary = self.getSummary(outputset)
        outputset.setObjComment(summary)
        outputs = {outputName: outputset}
        self._defineOutputs(**outputs)

        # Using a pointer to define the relations is more robust to scheduling
        # and id changes between the protocol run.db and the main project
        # database. The pointer defined below points to the outputset object
        self._defineSourceRelation(self.getInputMicrographsPointer(),
                                   Pointer(value=self, extended=outputName))
        self._store()


class ProtParticlePickingAuto(ProtParticlePicking):
    """ A derived class from ProtParticlePicking to differentiate those
    picking protocols that works in automatic way, i.e., that are not
    interactive and are also good candidates to be run in streaming. """

    def _insertAllSteps(self):
        self.initialIds = self._insertInitialSteps()
        self.micDict = OrderedDict()
        pwutils.makeFilePath(self._getAllDone())

        micDict, _ = self._loadInputList()
        pickMicIds = self._insertNewMicsSteps(micDict.values())

        self._insertFinalSteps(pickMicIds)

    def _insertInitialSteps(self):
        """ Override this function to insert some steps before the
        picking micrograph steps.
        Should return a list of ids of the initial steps. """
        return []

    def _insertFinalSteps(self, micSteps):
        """ Override this function to insert some steps after the
        picking micrograph steps.
        Receive the list of step ids of the picking steps. """
        self._insertFunctionStep('createOutputStep',
                                 prerequisites=micSteps, wait=True)

    def _getPickArgs(self):
        """ Should be implemented in sub-classes to define the argument
        list that should be passed to the picking step function.
        """
        return []

    def _insertPickMicrographStep(self, mic, prerequisites, *args):
        """ Basic method to insert a picking step for a given micrograph. """
        micStepId = self._insertFunctionStep('pickMicrographStep',
                                             mic.getMicName(), *args,
                                             prerequisites=prerequisites)

        return micStepId

    def pickMicrographStep(self, micName, *args):
        """ Step function that will be common for all picking protocols.
        It will take care of re-building the micrograph object from the micDict
        argument and perform any conversion if needed. Then, the function
        _pickMicrograph will be called, that should be implemented by each
        picking protocol.
        """
        mic = self.micDict[micName]
        micDoneFn = self._getMicDone(mic)
        micFn = mic.getFileName()

        if (self.isContinued() and os.path.exists(micDoneFn)):
            self.info("Skipping micrograph: %s, seems to be done" % micFn)
            return

        # Clean old finished files
        pwutils.cleanPath(micDoneFn)

        self.info("Picking micrograph: %s " % micFn)
        self._pickMicrograph(mic, *args)

        # Mark this mic as finished
        open(micDoneFn, 'w').close()

    def _pickMicrograph(self, mic, *args):
        """ This function should be implemented by subclasses in order
        to picking the given micrograph. """
        pass

    # --------------------------- UTILS functions ----------------------------

    # ------ Methods for Streaming picking --------------

    def _stepsCheck(self):
        # To allow streaming picking we need to detect:
        #   1) new micrographs ready to be picked
        #   2) new output coordinates that have been produced and add then
        #      to the output set.
        self._checkNewInput()
        self._checkNewOutput()

    def _insertNewMicsSteps(self, inputMics):
        """ Insert steps to process new mics (from streaming)
        Params:
            insertedDict: contains already processed mics
            inputMics: input mics set to be check
        """
        deps = []
        # For each mic insert the step to process it
        for mic in inputMics:
            micKey = mic.getMicName()
            if micKey not in self.micDict:
                stepId = self._insertPickMicrographStep(mic, self.initialIds,
                                                        *self._getPickArgs())
                deps.append(stepId)
                self.micDict[micKey] = mic
        return deps

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

    def _loadMics(self, micSet):
        return self._loadSet(micSet, SetOfMicrographs,
                        lambda mic: mic.getMicName())

    def _loadCTFs(self, ctfSet):
        return self._loadSet(ctfSet, SetOfCTF,
                        lambda ctf: ctf.getMicrograph().getMicName())

    def _loadInputList(self):
        """ Load the input set of micrographs that are ready to be picked. """
        return self._loadMics(self.getInputMicrographs())

    def _checkNewInput(self):
        # Check if there are new micrographs to process from the input set
        localFile = self.getInputMicrographs().getFileName()
        now = datetime.now()
        self.lastCheck = getattr(self, 'lastCheck', now)
        mTime = datetime.fromtimestamp(os.path.getmtime(localFile))
        self.debug('Last check: %s, modification: %s'
                  % (pwutils.prettyTime(self.lastCheck),
                     pwutils.prettyTime(mTime)))
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
            newDoneUpdated = self._updateOutputCoordSet(newDone, streamMode)
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

    def _micIsReady(self, mic):
        """ Function to check if a micrograph (although reported done)
        is ready to update the coordinates from it. An practical use of this
        function will be for protocols that need to wait for the CTF of that
        micrograph to be ready as well.
        """
        return True

    def readCoordsFromMics(self, outputDir, micDoneList, outputCoords):
        """ This method should be implemented in subclasses to read
        the coordinates from a given list of micrographs.
        """
        pass

    def _updateOutputCoordSet(self, micList, streamMode):
        micDoneList = [mic for mic in micList if self._micIsReady(mic)]

        # Do no proceed if there is not micrograph ready
        if not micDoneList:
            return []

        outputName = 'outputCoordinates'
        outputDir = self.getCoordsDir()
        outputCoords = getattr(self, outputName, None)

        # If there are not outputCoordinates yet, it means that is the first
        # time we are updating output coordinates, so we need to first create
        # the output set
        firstTime = outputCoords is None

        if firstTime:
            micSetPtr = self.getInputMicrographsPointer()
            outputCoords = self._createSetOfCoordinates(micSetPtr)
        else:
            outputCoords.enableAppend()

        self.readCoordsFromMics(outputDir, micDoneList, outputCoords)
        self.debug(" _updateOutputCoordSet Stream Mode: %s " % streamMode)
        self._updateOutputSet(outputName, outputCoords, streamMode)

        if firstTime:
            self._defineSourceRelation(self.getInputMicrographsPointer(),
                                       outputCoords)

        return micDoneList

    def _updateStreamState(self, streamMode):

        outputName = 'outputCoordinates'
        outputCoords = getattr(self, outputName, None)

        # If there are not outputCoordinates yet, it means that is the first
        # time we are updating output coordinates, so we need to first create
        # the output set
        firstTime = outputCoords is None

        if firstTime:
            micSetPtr = self.getInputMicrographsPointer()
            outputCoords = self._createSetOfCoordinates(micSetPtr)
        else:
            outputCoords.enableAppend()

        self.debug(" _updateStreamState Stream Mode: %s " % streamMode)
        self._updateOutputSet(outputName, outputCoords, streamMode)

    def _getMicDone(self, mic):
        return self._getExtraPath('DONE', 'mic_%06d.TXT' % mic.getObjId())

    def _isMicDone(self, mic):
        """ A mic is done if the marker file exists. """
        return os.path.exists(self._getMicDone(mic))

    def _getAllDone(self):
        return self._getExtraPath('DONE', 'all.TXT')

    def _readDoneList(self):
        """ Read from a text file the id's of the items that have been done. """
        doneFile = self._getAllDone()
        doneList = []
        # Check what items have been previously done
        if os.path.exists(doneFile):
            with open(doneFile) as f:
                doneList += [int(line.strip()) for line in f]

        return doneList

    def _writeDoneList(self, micList):
        """ Write to a text file the items that have been done. """
        doneFile = self._getAllDone()

        if not os.path.exists(doneFile):
            pwutils.makeFilePath(doneFile)

        with open(doneFile, 'a') as f:
            for mic in micList:
                f.write('%d\n' % mic.getObjId())

    def _getFirstJoinStepName(self):
        # This function will be used for streaming, to check which is
        # the first function that need to wait for all micrographs
        # to have completed, this can be overwritten in subclasses
        # (eg in Xmipp 'sortPSDStep')
        return 'createOutputStep'

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None

    def createOutputStep(self):
        # Not really required now
        #self._createOutput(self._getExtraPath())
        pass


# Micrograph type constants for particle extraction
SAME_AS_PICKING = 0
OTHER = 1

class ProtExtractParticles(ProtParticles):
    """ Base class for all extract-particles protocols.
     This class will take care of the streaming functionality and
     derived classes should mainly overwrite the '_extractMicrograph' function.
     """
    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
    
        form.addParam('inputCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates',
                      important=True,
                      label="Input coordinates",
                      help='Select the SetOfCoordinates ')
    
        # The name for the followig param is because historical reasons
        # now it should be named better 'micsSource' rather than
        # 'downsampleType', but this could make inconsistent previous executions
        # of this protocols, we will keep the name
        form.addParam('downsampleType', params.EnumParam,
                      choices=['same as picking', 'other'],
                      default=0, important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Micrographs source',
                      help='By default the particles will be extracted '
                           'from the micrographs used in the picking '
                           'step ( _same as picking_ option ). \n'
                           'If you select _other_ option, you must provide '
                           'a different set of micrographs to extract from. \n'
                           '*Note*: In the _other_ case, ensure that provided '
                           'micrographs and coordinates are related '
                           'by micName or by micId. Difference in pixel size '
                           'will be handled automatically.')
    
        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      condition='downsampleType != %s' % SAME_AS_PICKING,
                      important=True, label='Input micrographs',
                      help='Select the SetOfMicrographs from which to extract.')
    
        form.addParam('ctfRelations', params.RelationParam, allowsNull=True,
                      relationName=RELATION_CTF,
                      attributeName='getInputMicrographs',
                      label='CTF estimation',
                      help='Choose some CTF estimation related to input '
                           'micrographs. \n CTF estimation is needed if you '
                           'want to do phase flipping or you want to '
                           'associate CTF information to the particles.')
        
        self._definePreprocessParams(form)

    def _definePreprocessParams(self, form):
        """ Should be implemented in sub-classes to define some
        specific parameters """
        pass
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        # Let's load input data for the already existing micrographs
        # before the streaming
        self.debug(">>> _insertAllSteps ")
        pwutils.makeFilePath(self._getAllDone())
        self.micDict = OrderedDict()
        self.coordDict = {}

        micDict = self._loadInputList()

        self.initialIds = self._insertInitialSteps()

        pickMicIds = self._insertNewMicsSteps(micDict.values())

        self._insertFinalSteps(pickMicIds)

    def _insertInitialSteps(self):
        """ Override this function to insert some steps before the
        extract micrograph steps.
        Should return a list of ids of the initial steps. """
        return []

    def _insertNewMicsSteps(self, inputMics):
        """ Insert steps to process new mics (from streaming)
        Params:
            inputMics: input mics set to be check
        """
        deps = []
        #  TODO: We must check if the new inserted micrographs has associated
        #  coordinates. If not, we cannot extract particles. This is mandatory
        # if the inputMics are from "other" option.
    
        # For each mic insert the step to process it
        for mic in inputMics:
            micKey = mic.getMicName()
            if micKey not in self.micDict:
                args = self._getExtractArgs()
                stepId = self._insertExtractMicrographStep(mic, self.initialIds,
                                                           *args)
                deps.append(stepId)
                self.micDict[micKey] = mic
        return deps

    def _insertFinalSteps(self, micSteps):
        """ Override this function to insert some steps after the
        extraction of micrograph steps.
        Receive the list of step ids of the picking steps. """
        self._insertFunctionStep('createOutputStep',
                                 prerequisites=micSteps, wait=True)

    def _insertExtractMicrographStep(self, mic, prerequisites, *args):
        """ Basic method to insert a picking step for a given micrograph. """
        micStepId = self._insertFunctionStep('extractMicrographStep',
                                             mic.getMicName(), *args,
                                             prerequisites=prerequisites)
        return micStepId
    
    #--------------------------- STEPS functions -------------------------------
    def extractMicrographStep(self, micKey, *args):
        """ Step function that will be common for all extraction protocols.
        It will take an id and will grab the micrograph from a micDict map.
        The micrograph will be passed as input to the _extractMicrograph
        function.
        """
        # Retrieve the corresponding micrograph with this key and the
        # associated list of coordinates
        mic = self.micDict[micKey]
        coordList = self.coordDict[mic.getObjId()]
        self._convertCoordinates(mic, coordList)

        micDoneFn = self._getMicDone(mic)
        micFn = mic.getFileName()

        if (self.isContinued() and os.path.exists(micDoneFn)):
            self.info("Skipping micrograph: %s, seems to be done" % micFn)
            return

        # Clean old finished files
        pwutils.cleanPath(micDoneFn)

        self.info("Extracting micrograph: %s " % micFn)
        self._extractMicrograph(mic, *args)

        # Mark this mic as finished
        open(micDoneFn, 'w').close()

    def _extractMicrograph(self, mic, *args):
        """ This function should be implemented by subclasses in order
        to picking the given micrograph. """
        pass

    # --------------------------- UTILS functions ------------------------------
    def _convertCoordinates(self, mic, coordList):
        """ This function should be implemented by subclasses. """
        pass

    def _getExtractArgs(self):
        """ Should be implemented in sub-classes to define the argument
        list that should be passed to the extract step function.
        """
        return []

    def _micsOther(self):
        """ Return True if other micrographs are used for extract.
        Should be implemented in derived classes.
        """
        return False

    def _useCTF(self):
        """ Return True if a SetOfCTF is associated with the extracted
        particles.
         Should be implemented in derived classes.
        """
        return False
    
    # ------ Methods for Streaming extraction --------------
    def _isStreamClosed(self):
        return self.micsClosed and self.ctfsClosed and self.coordsClosed

    def _stepsCheck(self):
        # To allow streaming picking we need to detect:
        #   1) new micrographs ready to be picked
        #   2) new output coordinates that have been produced and add then
        #      to the output set.
        self._checkNewInput()
        self._checkNewOutput()

    def _loadInputList(self):
        """ Load the input set of micrographs and create a dictionary.
        The dictionary self.micDict, will contain all the micrographs
        from which particles have been extracted and also ones that
        are ready to be extracted.
        For creating this dictionary, we need to inspect:
        1) Input micrographs coming from coordinates, to know that a given
        micrographs has been picked.
        2) Micrographs to be extracted (in case it is different from 1)
        3) New computed CTF (in case it is associated with the particles)
        """
        def _loadSet(inputSet, SetClass, getKeyFunc):
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

        def _loadMics(micSet):
            return _loadSet(micSet, SetOfMicrographs,
                             lambda mic: mic.getMicName())

        def _loadCTFs(ctfSet):
            return _loadSet(ctfSet, SetOfCTF,
                             lambda ctf: ctf.getMicrograph().getMicName())

        # Load new micrographs coming from the coordinates
        self.debug("Loading Mics from Coords.")
        coordMics = self.inputCoordinates.get().getMicrographs()
        micDict, self.micsClosed = _loadMics(coordMics)
        # If we are extracting from other micrographs, then we will use
        # the other micrographs and filter those that
        if self._micsOther():
            self.debug("Loading other Mics.")
            oMicDict, oMicClosed = _loadMics(self.inputMicrographs.get())
            self.micsClosed = self.micsClosed and oMicClosed
            for micKey, mic in micDict.iteritems():
                if micKey in oMicDict:
                    oMic = oMicDict[micKey]
                    # Let's fix the id in case it does not correspond
                    # we want to have the id coming from the coordinates
                    # to match each coordinate to its micrograph
                    oMic.copyObjId(mic)
                    micDict[micKey] = oMic
                else:
                    del micDict[micKey]

        self.debug("Mics are closed? %s" % self.micsClosed)
        if self._useCTF():
            self.debug("Loading CTFs.")
            ctfDict, ctfClosed = _loadCTFs(self.ctfRelations.get())
            for micKey, mic in micDict.iteritems():
                if micKey in ctfDict:
                    mic.setCTF(ctfDict[micKey])
                else:
                    del micDict[micKey]
        # if not use CTF, self.ctfsClosed is True
        self.ctfsClosed = ctfClosed if self._useCTF() else True
        self.debug("CTFs are closed? %s" % self.ctfsClosed)

        self.debug("Loading Coords.")
        # Now load the coordinates for the newly detected micrographs. If
        # microgrpahs does not have coordinates, is not processed.
        micDict = self._loadInputCoords(micDict)

        return micDict

    def _loadInputCoords(self, micDict):
        """ Load coordinates from the input streaming.
        """
        coordsFn = self.getCoords().getFileName()
        self.debug("Loading input db: %s" % coordsFn)
        coordSet = SetOfCoordinates(filename=coordsFn)
        # FIXME: Temporary to avoid loadAllPropertiesFail
        coordSet._xmippMd = String()
        coordSet.loadAllProperties()

        for micKey, mic in micDict.iteritems():
            micId = mic.getObjId()
            coordList = []
            self.debug("Loading coords for mic: %s (%s)" % (micId, micKey))
            for coord in coordSet.iterItems(where='_micId=%s' % micId):
                # TODO: Check performance penalty of using this clone
                coordList.append(coord.clone())
            self.debug("   Coords found: %s" % len(coordList))

            if coordList:
                self.coordDict[micId] = coordList
            else:
                del micDict[micKey]
        self.coordsClosed = coordSet.isStreamClosed()
        coordSet.close()
        self.debug("Coords are closed? %s" % self.coordsClosed)
        self.debug("Closed db.")

        return micDict

    def _checkNewInput(self):
        self.debug(">>> _checkNewInput ")
    
        def _modificationTime():
            """ Check the last modification time of any of the three possible
             input files. """
            items = [self.inputCoordinates.get()]
        
            if self._micsOther():
                items.append(self.inputMicrographs.get())
            else:
                items.append(self.inputCoordinates.get().getMicrographs())
        
            if self._useCTF():
                items.append(self.ctfRelations.get())
        
            def _mTime(fn):
                return datetime.fromtimestamp(os.path.getmtime(fn))
        
            return max([_mTime(i.getFileName()) for i in items])
    
        mTime = _modificationTime()
        now = datetime.now()
        self.lastCheck = getattr(self, 'lastCheck', now)
        self.debug('Last check: %s, modification: %s'
                   % (pwutils.prettyTime(self.lastCheck),
                      pwutils.prettyTime(mTime)))
        # If the input micrographs.sqlite have not changed since our last check,
        # it does not make sense to check for new input data, but we must
        # check if sets are closed.
        self.debug("self.lastCheck > mTime %s , hasattr(self, 'micDict') %s"
                   % (self.lastCheck > mTime, hasattr(self, 'micDict')))
        if self.lastCheck > mTime and hasattr(self, 'micDict'):
            return None
        
        # Open input micrographs.sqlite and close it as soon as possible
        newMics = self._loadInputList()
        
        self.lastCheck = now
        outputStep = self._getFirstJoinStep()
    
        if newMics:
            fDeps = self._insertNewMicsSteps(newMics.values())
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()

    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return

        # Load previously done items (from text file)
        doneList = self._readDoneList()
        # Check for newly done items
        newDone = [m for m in self.micDict.values()
                   if m.getObjId() not in doneList and self._isMicDone(m)]

        # Update the file with the newly done mics
        # or exit from the function if no new done mics
        inputLen = len(self.micDict)
        self.debug('_checkNewOutput: ')
        self.debug('   input: %s, doneList: %s, newDone: %s'
                   % (inputLen, len(doneList), len(newDone)))

        firstTime = len(doneList) == 0
        allDone = len(doneList) + len(newDone)
        # We have finished when there is not more input mics (stream closed)
        # and the number of processed mics is equal to the number of inputs
        streamClosed = self._isStreamClosed()
        self.finished = streamClosed and allDone == inputLen
        self.debug(' is finished? %s ' % self.finished)
        self.debug(' is stream closed? %s ' % streamClosed)
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN
        if newDone:
            self._updateOutputPartSet(newDone, streamMode)
            self._writeDoneList(newDone)
        elif not self.finished:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        self.debug('   finished: %s ' % self.finished)
        self.debug('        self.streamClosed (%s) AND' % streamClosed)
        self.debug('        allDone (%s) == len(self.listOfMics (%s)'
                   % (allDone, inputLen))
        self.debug('   streamMode: %s' % streamMode)

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)

    def readCoordsFromMics(self, outputDir, micDoneList, outputCoords):
        """ This method should be implemented in subclasses to read
        the coordinates from a given list of micrographs.
        """
        pass

    def _updateOutputPartSet(self, micList, streamMode):
        outputName = 'outputParticles'
        outputParts = getattr(self, outputName, None)
        firstTime = True

        if outputParts is None:
            inputMics = self.getInputMicrographs()
            outputParts = self._createSetOfParticles()
            outputParts.setCoordinates(self.getCoords())
            if self.doFlip:
                outputParts.setIsPhaseFlipped(not inputMics.isPhaseFlipped())

            outputParts.setSamplingRate(self._getNewSampling())
            outputParts.setHasCTF(self._useCTF())
        else:
            firstTime = False
            outputParts.enableAppend()

        self.readPartsFromMics(micList, outputParts)
        self._updateOutputSet(outputName, outputParts, streamMode)

        if firstTime:
            # self._storeMethodsInfo(fnImages)
            self._defineSourceRelation(self.inputCoordinates, outputParts)
            if self._useCTF():
                self._defineSourceRelation(self.ctfRelations, outputParts)
            if self._micsOther():
                self._defineSourceRelation(self.inputMicrographs, outputParts)

    def _getMicDone(self, mic):
        return self._getExtraPath('DONE', 'mic_%06d.TXT' % mic.getObjId())

    def _isMicDone(self, mic):
        """ A mic is done if the marker file exists. """
        return os.path.exists(self._getMicDone(mic))

    def _getAllDone(self):
        return self._getExtraPath('DONE', 'all.TXT')

    def _readDoneList(self):
        """ Read from a text file the id's of the items that have been done. """
        doneFile = self._getAllDone()
        doneList = []
        # Check what items have been previously done
        if os.path.exists(doneFile):
            with open(doneFile) as f:
                doneList += [int(line.strip()) for line in f]

        return doneList

    def _writeDoneList(self, micList):
        """ Write to a text file the items that have been done. """
        doneFile = self._getAllDone()

        if not os.path.exists(doneFile):
            pwutils.makeFilePath(doneFile)

        with open(doneFile, 'a') as f:
            for mic in micList:
                f.write('%d\n' % mic.getObjId())

    def _getFirstJoinStepName(self):
        # This function will be used for streaming, to check which is
        # the first function that need to wait for all micrographs
        # to have completed, this can be overwritten in subclasses
        # (eg in Xmipp 'sortPSDStep')
        return 'createOutputStep'

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None

    def createOutputStep(self):
        pass # Nothing to do now
        #self._createOutput(self._getExtraPath())


class ProtExtractParticlesPair(ProtParticles):
    """ Base class for all extract-particles pairs protocols. Until now,
    this protcols is not in streaming mode.
     """
