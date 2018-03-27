# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
import time
from datetime import datetime
from collections import OrderedDict

from pyworkflow.object import Set, Pointer
import pyworkflow.protocol.params as params
from pyworkflow.protocol import STATUS_NEW
from pyworkflow.em.data import (EMObject, SetOfCoordinates, SetOfMicrographs,
                                SetOfCTF)
import pyworkflow.utils as pwutils
from pyworkflow.utils.properties import Message

from protocol_particles import ProtParticles


class ProtParticlePicking(ProtParticles):
    OUTPUT_PREFIX = 'outputCoordinates'

    def _defineParams(self, form):

        form.addSection(label='Input')
        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      label=Message.LABEL_INPUT_MIC, important=True,
                      help='Select the SetOfMicrographs to be used during '
                           'picking.')

    # -------------------------- INFO functions -------------------------------
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

        micDict, self.streamClosed = self._loadInputList()
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

        if self.isContinued() and os.path.exists(micDoneFn):
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

    # Group of functions to pick several micrographs if the batch size is
    # defined. In some programs it might be more efficient to pick many
    # at once and not one by one

    def _insertPickMicrographListStep(self, micList, prerequisites, *args):
        """ Basic method to insert a picking step for a given micrograph. """
        micNameList = [mic.getMicName() for mic in micList]
        micStepId = self._insertFunctionStep('pickMicrographListStep',
                                             micNameList, *args,
                                             prerequisites=prerequisites)
        return micStepId

    def pickMicrographListStep(self, micNameList, *args):
        micList = []

        for micName in micNameList:
            mic = self.micDict[micName]
            micDoneFn = self._getMicDone(mic)
            micFn = mic.getFileName()
            if self.isContinued() and os.path.exists(micDoneFn):
                self.info("Skipping micrograph: %s, seems to be done" % micFn)

            else:
                # Clean old finished files
                pwutils.cleanPath(micDoneFn)
                self.info("Picking micrograph: %s " % micFn)
                micList.append(mic)

        self._pickMicrographList(micList, *args)

        for mic in micList:
            # Mark this mic as finished
            open(self._getMicDone(mic), 'w').close()

    def _pickMicrographList(self, micList, *args):
        """ This function can be implemented by subclasses if it is a more
        effcient way to pick many micrographs at once.
         Default implementation will just call the _pickMicrograph
        """
        for mic in micList:
            self._pickMicrograph(mic, *args)

    # --------------------------- UTILS functions ----------------------------

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
        return self._insertNewMics(inputMics,
                                   lambda mic: mic.getMicName(),
                                   self._insertPickMicrographStep,
                                   self._insertPickMicrographListStep,
                                   *self._getPickArgs())

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

            # Maybe it would be good idea to take a snap to avoid
            # so much IO if this protocol does not have much to do now
            if allDone == nMics:
                self._streamingSleepOnWait()

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

        self.info("Reading coordinates from mics: %s" % ','.join([mic.strId() for mic in micList]))
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


class ProtPickingDifference(ProtParticlePicking):
    """
    Protocol to compute the difference between a reference SetOfPartices and
    a another set (usually a negative reference).

    The output will be a SetOfCoordinates with the particles in the reference
    input that are not close to coordinates in the negative input.

    """
    _label = 'picking difference'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates',
                      label="Input coordinates",
                      help='Select the reference set of coordinates.')
        form.addParam('negativeCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates',
                      label="Negative coordinates",
                      help='Negative coordinates that will help to exclude '
                           'coordinates from the input set. ')
        form.addParam('differenceRadius', params.IntParam, default=10,
                      label="Radius (px)",
                      help="Distance radius to consider that two coordinates "
                           "close enough. If a particle in the input reference"
                           "set have any close particle in the negative set, "
                           "it will not be included in the output set. ")

        form.addParallelSection(threads=0, mpi=0)

    # ------------------------- INSERT steps functions -------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep("createOutputStep",
                                 self.inputCoordinates.getObjId(),
                                 self.negativeCoordinates.getObjId(),
                                 self.differenceRadius.get())

    def getInputMicrographs(self):
        return self.inputCoordinates[0].get().getMicrographs()

    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        return []

    def createOutputStep(self, inputId, negId, radius):
        inputCoords = self.inputCoordinates.get()
        negCoords = self.negativeCoordinates.get()
        inputMics = inputCoords.getMicrographs()
        outputCoords = self._createSetOfCoordinates(inputMics)
        outputCoords.setBoxSize(inputCoords.getBoxSize())

        r2 = radius * radius  # to avoid computing sqrt when comparing distances

        def _discard(coord, negPos):
            """ Find if there is a close negative particle to this one. """
            x, y = coord.getPosition()

            for nx, ny in negPos:
                if abs((x - nx) * (y - ny)) < r2:
                    return True

            return False  # Far enough from all negative coordinates

        # Read all consensus particles
        for mic in inputMics:
            negPos = [(c.getX(), c.getY())
                      for c in negCoords.iterCoordinates(mic)]
            for coord in inputCoords.iterCoordinates(mic):
                if not _discard(coord, negPos):
                    outputCoords.append(coord)

        # Set output
        self._defineOutputs(outputCoordinates=outputCoords)
        self._defineTransformRelation(self.inputCoordinates,
                                      self.outputCoordinates)
        self._defineSourceRelation(self.negativeCoordinates,
                                   self.outputCoordinates)