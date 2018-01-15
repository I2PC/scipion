# **************************************************************************
# *
# * Authors:     Tomas Majtner (tmajtner@cnb.csic.es)
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
import pyworkflow.em as em
import pyworkflow.protocol.constants as cons
from pyworkflow.em.data import SetOfParticles
from pyworkflow.em.protocol import ProtProcessParticles
from pyworkflow.object import Set
from pyworkflow.protocol.params import (EnumParam, IntParam, Positive, Range,
                                        LEVEL_ADVANCED, FloatParam,
                                        BooleanParam)
from pyworkflow.em.packages.xmipp3.convert import readSetOfParticles, \
    writeSetOfParticles, setXmippAttributes


class XmippProtTriggerData(ProtProcessParticles):
    """ Waits until certain number of particles is prepared and then
    send them to output. """

    _label = 'trigger data'

    # --------------------------- DEFINE param functions ----------------------
    def _defineProcessParams(self, form):

        form.addParam('outputSize', IntParam, default=10000,
                      label='Output size',
                      help='How many particles need to be on input to '
                           'create output set.')
        form.addParam('checkInterval', IntParam, default=60,
                      label="Check for new objects each (sec)",
                      help="Check for new objects each checkInterval seconds")
        form.addParam('allParticles', BooleanParam, default=False,
                      label='Send all particles to output?',
                      help='If NO is selected, only subset of "Output size" '
                           'particles will be send to output.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self.SetOfParticles = [m.clone() for m in self.inputParticles.get()]
        partsSteps = self._checkNewPartsSteps()
        self._insertFunctionStep('createOutputStep',
                                 prerequisites=partsSteps, wait=True)

    def createOutputStep(self):
        pass

    def _getFirstJoinStepName(self):
        # This function will be used for streaming, to check which is
        # the first function that need to wait for all micrographs
        # to have completed, this can be overriden in subclasses
        # (e.g., in Xmipp 'sortPSDStep')
        return 'createOutputStep'

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None

    def _insertNewPartsSteps(self, insertedDict, inputParts):
        deps = []
        writeSetOfParticles([p.clone() for p in inputParts],
                            self._getExtraPath("allDone.xmd"),
                            alignType=em.ALIGN_NONE)
        stepId = self._insertFunctionStep('checkForOutput',
                                          inputParts,
                                          prerequisites=[])
        deps.append(stepId)
        for p in inputParts:
            if int(p.getObjId()) not in insertedDict:
                insertedDict[p.getObjId()] = stepId
        return deps

    def _stepsCheck(self):
        # Input particles set can be loaded or None when checked for new inputs
        # If None, we load it
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        # Check if there are new particles to process from the input set
        partsFile = self.inputParticles.get().getFileName()
        partsSet = SetOfParticles(filename=partsFile)
        partsSet.loadAllProperties()
        self.SetOfParticles = [m.clone() for m in partsSet]
        self.streamClosed = partsSet.isStreamClosed()
        partsSet.close()
        partsSet = self._createSetOfParticles()
        readSetOfParticles(self._getExtraPath("allDone.xmd"), partsSet)
        newParts = any(m.getObjId() not in partsSet
                       for m in self.SetOfParticles)
        outputStep = self._getFirstJoinStep()
        if newParts:
            fDeps = self._insertNewPartsSteps(self.insertedDict,
                                              self.SetOfParticles)
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()

    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return
        # Load previously done items (from text file)
        doneList = self._readDoneList()
        # Check for newly done items
        partsSet = self._createSetOfParticles()
        readSetOfParticles(self._getExtraPath("allDone.xmd"), partsSet)
        newDone = [m.clone() for m in self.SetOfParticles
                   if m.getObjId() not in doneList]
        self.finished = (self.streamClosed
                         and (len(doneList) == len(partsSet))) \
                        or (len(doneList) > self.outputSize
                            and self.allParticles == False)
        if newDone:
            self._writeDoneList(newDone)
        elif not self.finished:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return
        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(cons.STATUS_NEW)

    def _loadOutputSet(self, SetClass, baseName, newParts):
        setFile = self._getPath(baseName)
        if os.path.exists(setFile):
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        inputs = self.inputParticles.get()
        outputSet.copyInfo(inputs)
        outputSet.copyItems(newParts)
        return outputSet

    def _updateOutputSet(self, outputName, outputSet, state=Set.STREAM_OPEN):
        outputSet.setStreamState(state)
        if self.hasAttribute(outputName):
            outputSet.write()  # Write to commit changes
            outputAttr = getattr(self, outputName)
            # Copy the properties to the object contained in the protocol
            outputAttr.copy(outputSet, copyId=False)
            # Persist changes
            self._store(outputAttr)
        else:
            self._defineOutputs(**{outputName: outputSet})
            self._store(outputSet)

        # Close set databaset to avoid locking it
        outputSet.close()

    # --------------------------- STEPS functions -----------------------------
    def checkForOutput(self, newParts):

        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN
        outSet = self._loadOutputSet(SetOfParticles,
                                     'outputParticles.sqlite',
                                     newParts)
        self._updateOutputSet('outputParticles', outSet, streamMode)
        if (outSet.getSize() >= self.outputSize):
            self.finished = True

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Not enough particles for output yet.")
        else:
            summary.append("Particles were send to output.")
        return summary

    def _validate(self):
        pass

    # --------------------------- UTILS functions -----------------------------
    def _readDoneList(self):
        """ Read from a file the id's of the items that have been done. """
        doneFile = self._getAllDone()
        doneList = []
        # Check what items have been previously done
        if os.path.exists(doneFile):
            with open(doneFile) as f:
                doneList += [int(line.strip()) for line in f]
        return doneList

    def _getAllDone(self):
        return self._getExtraPath('DONE_all.TXT')

    def _writeDoneList(self, partList):
        """ Write to a file the items that have been done. """
        with open(self._getAllDone(), 'a') as f:
            for part in partList:
                f.write('%d\n' % part.getObjId())