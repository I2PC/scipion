# *****************************************************************************
# *
# * Authors:     Tomas Majtner         tmajtner@cnb.csic.es (2017)
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
# *****************************************************************************

import os
import pyworkflow.em as em
import pyworkflow.em.metadata as md
import pyworkflow.protocol.constants as cons
import pyworkflow.protocol.params as param
from datetime import datetime
from pyworkflow.em.protocol import ProtClassify2D
from pyworkflow.em.data import SetOfParticles
from pyworkflow.object import Set
from pyworkflow.utils.properties import Message
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles, \
    readSetOfParticles, setXmippAttributes


class XmippProtEliminateEmptyParticles(ProtClassify2D):
    """ Takes a set of particles and using statistical methods (variance of
    variances of sub-parts of input image) eliminates those samples, where
    there is no object/particle (only noise is presented there). Threshold
    parameter can be used for fine-tuning the algorithm for type of data. """

    _label = 'eliminate empty particles'

    def __init__(self, **args):
        ProtClassify2D.__init__(self, **args)
        self.stepsExecutionMode = em.STEPS_PARALLEL

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('inputParticles', param.PointerParam,
                      label="Input images",
                      important=True, pointerClass='SetOfParticles',
                      help='Select the input images to be classified.')
        form.addParam('threshold', param.FloatParam, default=1.5,
                      label='Threshold used in elimination:',
                      help='Higher threshold => '
                           'more particles will be eliminated. '
                           'Set to -1 for no elimination.')
        form.addParam('addFeatures', param.BooleanParam, default=False,
                      label='Add features', expertLevel=param.LEVEL_ADVANCED,
                      help='Add features used for the ranking to each '
                           'one of the input particles')
        form.addParam('useDenoising', param.BooleanParam, default=False,
                      label='Turning on denoising',
                      expertLevel=param.LEVEL_ADVANCED,
                      help='Option for turning on denoising method '
                           'while computing emptiness feature')
        form.addParam('denoising', param.IntParam, default=50,
                      expertLevel=param.LEVEL_ADVANCED,
                      condition='useDenoising',
                      label='Denoising factor:',
                      help='Factor to be used during denoising operation.'
                           'Higher value applies stronger denoising,'
                           'could be more precise by also very slow.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self.outputSize = 0
        self.check = None
        self.fnOutputMd = self._getExtraPath("output.xmd")
        self.fnElimMd = self._getExtraPath("eliminated.xmd")
        self.fnOutMdTmp = self._getExtraPath("outTemp.xmd")
        self.fnElimMdTmp = self._getExtraPath("elimTemp.xmd")
        checkStep = self._insertNewPartsSteps()
        self._insertFunctionStep('createOutputStep',
                                 prerequisites=checkStep, wait=True)

    def createOutputStep(self):
        pass

    def _getFirstJoinStepName(self):
        # This function will be used for streaming, to check which is
        # the first function that need to wait for all particles
        # to have completed, this can be overriden in subclasses
        # (e.g., in Xmipp 'sortPSDStep')
        return 'createOutputStep'

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None

    def _insertNewPartsSteps(self):
        deps = []
        stepId = self._insertFunctionStep('eliminationStep',
                                          self._getExtraPath("input.xmd"),
                                          prerequisites=[])
        deps.append(stepId)
        return deps

    def eliminationStep(self, fnInputMd):
        partsFile = self.inputParticles.get().getFileName()
        self.partsSet = SetOfParticles(filename=partsFile)
        self.partsSet.loadAllProperties()
        self.streamClosed = self.partsSet.isStreamClosed()
        if self.check == None:
            writeSetOfParticles(self.partsSet, fnInputMd,
                                alignType=em.ALIGN_NONE, orderBy='creation')
        else:
            writeSetOfParticles(self.partsSet, fnInputMd,
                                alignType=em.ALIGN_NONE, orderBy='creation',
                                where='creation>"' + str(self.check) + '"')
        for p in self.partsSet.iterItems(orderBy='creation', direction='DESC'):
            self.check = p.getObjCreation()
            break
        self.partsSet.close()
        args = "-i %s -o %s -e %s -t %f" % (
        fnInputMd, self.fnOutputMd, self.fnElimMd, self.threshold.get())
        if self.addFeatures:
            args+=" --addFeatures"
        if self.useDenoising:
            args += " --useDenoising -d %d" % self.denoising.get()
        self.runJob("xmipp_image_eliminate_empty_particles", args)
        os.remove(fnInputMd)
        streamMode = Set.STREAM_CLOSED \
            if getattr(self, 'finished', False) else Set.STREAM_OPEN
        if os.path.exists(self.fnOutMdTmp):
            outSet = self._loadOutputSet(SetOfParticles,
                                         'outputParticles.sqlite',
                                         self.fnOutMdTmp)
            self._updateOutputSet('outputParticles', outSet, streamMode)
        if os.path.exists(self.fnElimMdTmp):
            outSet = self._loadOutputSet(SetOfParticles,
                                         'eliminatedParticles.sqlite',
                                         self.fnElimMdTmp)
            self._updateOutputSet('eliminatedParticles', outSet, streamMode)

    def _stepsCheck(self):
        # Input particles set can be loaded or None when checked for new inputs
        # If None, we load it
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        # Check if there are new particles to process from the input set
        partsFile = self.inputParticles.get().getFileName()
        self.lastCheck = getattr(self, 'lastCheck', datetime.now())
        mTime = datetime.fromtimestamp(os.path.getmtime(partsFile))
        # If the input movies.sqlite have not changed since our last check,
        # it does not make sense to check for new input data
        if self.lastCheck > mTime:
            return None
        outputStep = self._getFirstJoinStep()
        fDeps = self._insertNewPartsSteps()
        if outputStep is not None:
            outputStep.addPrerequisites(*fDeps)
        self.updateSteps()

    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return
        self.finished = self.streamClosed and \
                        self.outputSize == len(self.partsSet)
        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(cons.STATUS_NEW)

    def _loadOutputSet(self, SetClass, baseName, fnMd):
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
        partsSet = self._createSetOfParticles()
        readSetOfParticles(fnMd, partsSet)
        self.outputSize = self.outputSize + len(partsSet)
        outputSet.copyItems(partsSet,
                          updateItemCallback=self._updateParticle,
                          itemDataIterator=
                            md.iterRows(fnMd, sortByLabel=md.MDL_ITEM_ID))
        os.remove(fnMd)
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

    # --------------------------- UTILS functions -----------------------------
    def _updateParticle(self, item, row):
        setXmippAttributes(item, row, md.MDL_SCORE_BY_EMPTINESS)
        if row.getValue(md.MDL_ENABLED) <= 0:
            item._appendItem = False
        else:
            item._appendItem = True
