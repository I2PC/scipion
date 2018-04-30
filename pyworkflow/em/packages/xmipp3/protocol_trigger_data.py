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
import time
from datetime import datetime
import pyworkflow.protocol.constants as cons
from pyworkflow.em.data import SetOfParticles
from pyworkflow.em.protocol import ProtProcessParticles
from pyworkflow.object import Set
from pyworkflow.protocol.params import (BooleanParam, IntParam)

class XmippProtTriggerData(ProtProcessParticles):
    """ Waits until certain number of particles is prepared and then
    send them to output. """

    _label = 'trigger data'

    # --------------------------- DEFINE param functions ----------------------
    def _defineProcessParams(self, form):

        form.addParam('outputSize', IntParam, default=10000,
                      label='Minimum output size',
                      help='How many particles need to be on input to '
                           'create output set.')
        form.addParam('checkInterval', IntParam, default=60,
                      label="Check for new objects each (sec)",
                      help="Check for new objects each checkInterval seconds")
        form.addParam('allParticles', BooleanParam, default=False,
                      label='Send all particles to output?',
                      help='If NO is selected, only subset of "Output size" '
                           'particles will be send to output.')
        form.addParam('delay', IntParam, default=1,
                      label="Delay (sec)", expertLevel=cons.LEVEL_ADVANCED,
                      help="Delay in seconds before checking new output")

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self.finished = False
        self.particles = []
        partsSteps = self.delayStep()
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

    def delayStep(self):
        time.sleep(self.delay)

    def _stepsCheck(self):
        # Input particles set can be loaded or None when checked for new inputs
        # If None, we load it
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        partsFile = self.inputParticles.get().getFileName()
        self.lastCheck = getattr(self, 'lastCheck', datetime.now())
        mTime = datetime.fromtimestamp(os.path.getmtime(partsFile))
        # If the input movies.sqlite have not changed since our last check,
        # it does not make sense to check for new input data
        if self.lastCheck > mTime and hasattr(self, 'newParticles'):
            return None

        # Open input movies.sqlite and close it as soon as possible
        self.partsSet = SetOfParticles(filename=partsFile)
        self.partsSet.loadAllProperties()
        if len(self.particles) > 0:
            self.newParticles = [m.clone() for m in self.partsSet.iterItems(
                orderBy='creation',
                where='creation>"' + str(self.check) + '"')]
        else:
            self.newParticles = [m.clone() for m in self.partsSet]

        self.particles = self.particles + self.newParticles
        if len(self.newParticles) > 0:
            for p in self.partsSet.iterItems(orderBy='creation',
                                             direction='DESC'):
                self.check = p.getObjCreation()
                break

        self.lastCheck = datetime.now()
        self.streamClosed = self.partsSet.isStreamClosed()
        self.partsSet.close()
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if len(self.particles) >= self.outputSize:
            if self.allParticles:
                if not os.path.exists(self._getPath('particles.sqlite')):
                    imageSet = self._loadOutputSet(SetOfParticles,
                                                   'particles.sqlite',
                                                   self.particles)
                else:
                    imageSet = self._loadOutputSet(SetOfParticles,
                                                   'particles.sqlite',
                                                   self.newParticles)
                self._updateOutputSet('outputParticles', imageSet, streamMode)

            elif not os.path.exists(self._getPath('particles.sqlite')):
                imageSet = self._loadOutputSet(SetOfParticles,
                                               'particles.sqlite',
                                               self.particles)
                self._updateOutputSet('outputParticles', imageSet, streamMode)

    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return
        self.finished = (not self.allParticles and
                         len(self.particles) > self.outputSize) or \
                        (self.streamClosed and self.allParticles
                          and len(self.particles) == len(self.partsSet))
        outputStep = self._getFirstJoinStep()
        deps = []
        if self.finished:  # Unlock createOutputStep if finished all jobs
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(cons.STATUS_NEW)
        else:
            delayId = self._insertFunctionStep('delayStep', prerequisites=[])
            deps.append(delayId)

        if outputStep is not None:
            outputStep.addPrerequisites(*deps)
        self.updateSteps()

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