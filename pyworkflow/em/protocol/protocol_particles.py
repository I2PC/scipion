# **************************************************************************
# *
# * Authors:     Airen Zaldivar Peraza (azaldivar@cnb.csic.es)
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

import os
from datetime import datetime

from pyworkflow.object import Set
from pyworkflow.protocol.params import PointerParam
from pyworkflow.protocol import STATUS_NEW
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.data import (EMObject, SetOfCoordinates, Micrograph,
                                SetOfMicrographs)
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
        
        form.addParam('inputParticles', PointerParam,
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


class ProtExtractParticles(ProtParticles):
    pass


class ProtParticlePicking(ProtParticles):
    OUTPUT_PREFIX = 'outputCoordinates'

    def _defineParams(self, form):

        form.addSection(label='Input')
        form.addParam('inputMicrographs', PointerParam,
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
        self._defineSourceRelation(self.getInputMicrographsPointer(), outputset)
        self._store()


class ProtParticlePickingAuto(ProtParticlePicking):
    """ A derived class from ProtParticlePicking to differentiate those
    picking protocols that works in automatic way, i.e., that are not
    interactive and are also good candidates to be run in streaming. """

    def _insertAllSteps(self):
        self.initialIds = self._insertInitialSteps()
        self.insertedDict = {}

        pickMicIds = self._insertNewMicsSteps(self.insertedDict,
                                              self.getInputMicrographs())

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
        micDict = mic.getObjDict(includeBasic=True)
        micStepId = self._insertFunctionStep('pickMicrographStep',
                                             micDict, *args,
                                             prerequisites=prerequisites)

        return micStepId

    def pickMicrographStep(self, micDict, *args):
        """ Step function that will be common for all picking protocols.
        It will take care of re-building the micrograph object from the micDict
        argument and perform any conversion if needed. Then, the function
        _pickMicrograph will be called, that should be implemented by each
        picking protocol.
        """
        mic = Micrograph()
        mic.setAttributesFromDict(micDict, setBasic=True, ignoreMissing=True)

        # Clean old finished files
        pwutils.cleanPath(self._getMicDone(mic))

        self.info("Picking micrograph: %s " % mic.getFileName())
        self._pickMicrograph(mic, *args)

        # Mark this mic as finished
        open(self._getMicDone(mic), 'w').close()

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

    def _insertNewMicsSteps(self, insertedDict, inputMics):
        """ Insert steps to process new mics (from streaming)
        Params:
            insertedDict: contains already processed mics
            inputMics: input mics set to be check
        """
        deps = []
        # For each mic insert the step to process it
        for mic in inputMics:
            if mic.getObjId() not in insertedDict:
                stepId = self._insertPickMicrographStep(mic, self.initialIds,
                                                        *self._getPickArgs())
                deps.append(stepId)
                insertedDict[mic.getObjId()] = stepId
        return deps

    def _loadInputList(self):
        """ Load the input set of micrographs and create a list. """
        micsFile = self.getInputMicrographs().getFileName()
        self.debug("Loading input db: %s" % micsFile)
        micSet = SetOfMicrographs(filename=micsFile)
        micSet.loadAllProperties()
        self.listOfMics = [m.clone() for m in micSet]
        self.streamClosed = micSet.isStreamClosed()
        micSet.close()
        self.debug("Closed db.")

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
        self._loadInputList()
        newMics = any(m.getObjId() not in self.insertedDict
                        for m in self.listOfMics)
        outputStep = self._getFirstJoinStep()

        if newMics:
            fDeps = self._insertNewMicsSteps(self.insertedDict, self.listOfMics)
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()

    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return

        # Load previously done items (from text file)
        doneList = self._readDoneList()
        # Check for newly done items
        newDone = [m for m in self.listOfMics
                   if m.getObjId() not in doneList and self._isMicDone(m)]

        # Update the file with the newly done mics
        # or exit from the function if no new done mics
        self.debug('_checkNewOutput: ')
        self.debug('   listOfMics: %s, doneList: %s, newDone: %s'
                   % (len(self.listOfMics), len(doneList), len(newDone)))

        firstTime = len(doneList) == 0
        allDone = len(doneList) + len(newDone)
        # We have finished when there is not more input mics (stream closed)
        # and the number of processed mics is equal to the number of inputs
        self.finished = self.streamClosed and allDone == len(self.listOfMics)
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if newDone:
            self._writeDoneList(newDone)
        elif not self.finished:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        self.debug('   finished: %s ' % self.finished)
        self.debug('        self.streamClosed (%s) AND' % self.streamClosed)
        self.debug('        allDone (%s) == len(self.listOfMics (%s)'
                   % (allDone, len(self.listOfMics)))
        self.debug('   streamMode: %s' % streamMode)

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)

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