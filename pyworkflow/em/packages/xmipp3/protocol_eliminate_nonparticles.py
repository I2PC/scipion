# ******************************************************************************
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
# ******************************************************************************

import os
from os.path import basename
import pyworkflow.em as em
import pyworkflow.em.metadata as md
import pyworkflow.protocol.constants as cons
import pyworkflow.protocol.params as param
from pyworkflow.em.protocol import ProtClassify2D, SetOfClasses2D
from pyworkflow.em.data import SetOfParticles
from pyworkflow.object import Set
from pyworkflow.utils import getExt, replaceExt
from pyworkflow.utils.path import cleanPath
from pyworkflow.utils.properties import Message
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles, xmippToLocation


class XmippProtEliminateNonParticles(ProtClassify2D):
    """ Takes a set of particles and using statistical methods (variance of
    variances of sub-parts of input image) eliminates those samples, where
    is no object/particle (only noise is presented there). Threshold parameter
    can be used for fine-tuning the algorithm for type of data. """

    _label = 'eliminate non-particles'

    def __init__(self, **args):
        ProtClassify2D.__init__(self, **args)
        self.stepsExecutionMode = em.STEPS_PARALLEL

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('inputParticles', param.PointerParam,
                      label="Input images",
                      important=True, pointerClass='SetOfParticles',
                      help='Select the input images to be classified.')
        form.addParam('threshold', param.IntParam, default=1,
                      label='Threshold used in elimination:',
                      help='Higher threshold = more particles are eliminated.')

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        writeSetOfParticles(self.inputParticles.get(),
                            self._getExtraPath("images.xmd"),
                            alignType=em.ALIGN_NONE)
        # Convert input images if necessary
        self.insertedDict = {}
        self._insertNewPartsSteps(self.insertedDict, self.inputParticles.get())
        self._insertFunctionStep('createOutputStep', wait=True)

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

    def _insertNewPartsSteps(self, insertedDict, inputParts):
        deps = []
        stepId = self._insertStepsForParticles(
            self._getExtraPath("images.xmd"), self._getExtraPath("output.xmd"))
        deps.append(stepId)
        # stepId = self._insertStepsForParticles(self._getExtraPath("images2.xmd"), self._getExtraPath("output.xmd"))
        # deps.append(stepId)
        for part in inputParts:
            if part.getObjId() not in insertedDict:
                insertedDict[part.getObjId()] = stepId
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
        newParts = any(m.getObjId() not in self.insertedDict
                       for m in self.inputParticles.get())
        outputStep = self._getFirstJoinStep()
        if newParts:
            fDeps = self._insertNewPartsSteps(self.insertedDict,
                                              self.inputParticles.get())
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()

    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return
        # Load previously done items (from text file)
        doneList = self._readDoneList()
        # Check for newly done items
        newDone = [m.clone() for m in self.SetOfParticles
                   if int(m.getObjId()) not in doneList]
        # We have finished when there is not more input particles (stream closed)
        # and the number of processed particles is equal to the number of inputs
        self.finished = self.streamClosed and (len(doneList) + len(
            newDone)) == len(self.SetOfParticles)
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if newDone:
            self._writeDoneList(newDone)
        elif not self.finished:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        outSet = self._loadOutputSet(SetOfClasses2D, 'classes2D.sqlite')

        # for part in newDone:
        #     partOut = em.data.Particle()
        #     partOut.setObjId(part.getObjId())
        #     partOut.setFileName(self._getOutputParticle(part))
        #     outSet.append(partOut)
        #
        self._updateOutputSet('outputParticles', outSet, streamMode)

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(cons.STATUS_NEW)

    def _loadOutputSet(self, SetClass, baseName):
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
            inputs = self.inputParticles.get()
            # Here the defineOutputs function will call the write() method
            classes2DSet = self._createSetOfClasses2D(inputs)
            self._fillClassesFromLevel(classes2DSet)

            outputSet = {'outputClasses': classes2DSet}
            self._defineOutputs(**outputSet)
            self._defineSourceRelation(self.inputParticles, classes2DSet)

    def _insertStepsForParticles(self, inputPart, outputPart):
        classifyStepId = self._insertFunctionStep('fastClassifyStep',
                                                  inputPart,
                                                  outputPart)
        return classifyStepId

    def fastClassifyStep(self, fnInputMd, fnOutputMd):
        iteration = 0
        args = "-i %s -o %s -k %d" % (
        fnInputMd, fnOutputMd, self.threshold.get())
        self.runJob("xmipp_classify_fast_2d", args)

        cleanPath(self._getExtraPath("level_00"))
        blocks = md.getBlocksInMetaDataFile(fnOutputMd)
        fnDir = self._getExtraPath()
        # Gather all images in block
        for block in blocks:
            if block.startswith('class0'):
                args = "-i %s@%s --iter 5 --distance correlation --classicalMultiref --nref 1 --odir %s --oroot %s" % (
                block, fnOutputMd, fnDir, block)
                if iteration == 0:
                    args += " --nref0 1"
                else:
                    args += " --ref0 %s" % self._getExtraPath(
                        "level_00/%s_classes.stk" % block)
                self.runJob("xmipp_classify_CL2D", args,
                            numberOfMpi=max(2, self.numberOfMpi.get()))
                cleanPath(
                    self._getExtraPath("level_00/%s_classes.xmd" % block))

    # --------------------------- UTILS functions -----------------------------
    def _defArgsClassify(self, inputPart, outputPart):
        # Prepare arguments to call program: xmipp_classify_CL2D
        args = "-i %s -o %s -k %d" % (
        inputPart, outputPart, self.numberOfClasses.get())
        if False:
            args += ' -c %(nref0)s'
        return args

    def _readDoneList(self):
        """ Read from a text file the id's of the items that have been done. """
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
        """ Write to a text file the items that have been done. """
        with open(self._getAllDone(), 'a') as f:
            for part in partList:
                f.write('%d\n' % part.getObjId())

    def _isPartDone(self, part):
        """ A particle is done if the marker file exists. """
        return os.path.exists(self._getPartDone(part))

    def _getPartDone(self, part):
        fn = str(part.getObjId()) + "_" + part.getFileName()
        extFn = getExt(fn)
        if extFn != ".stk":
            fn = replaceExt(fn, "stk")
        return self._getExtraPath('%s' % basename(fn))

    # --------------------------- UTILS functions -------------------------------

    def _updateParticle(self, item, row):
        item.setClassId(row.getValue(md.MDL_REF))

    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            index, fn, _ = self._classesInfo[classId]
            item.setAlignment2D()
            rep = item.getRepresentative()
            rep.setLocation(index, fn)
            rep.setSamplingRate(self.inputParticles.get().getSamplingRate())

    def _loadClassesInfo(self, filename, blockId):
        """ Read some information about the produced 2D classes
        from the metadata file.
        """
        self._classesInfo = {}  # store classes info, indexed by class id
        mdClasses = md.MetaData(filename)
        for classNumber, row in enumerate(md.iterRows(mdClasses)):
            index, fn = xmippToLocation(row.getValue(md.MDL_IMAGE))
            self._classesInfo[blockId] = (index, fn, row.clone())

    def _fillClassesFromLevel(self, clsSet):
        """ Create the SetOfClasses2D from a given iteration. """
        blocks = md.getBlocksInMetaDataFile(self._getExtraPath('output.xmd'))
        for blockId, block in enumerate(blocks):
            if block.startswith('class0'):
                self._loadClassesInfo(
                    self._getExtraPath("level_00/%s_classes.stk" % block),
                    blockId)
                xmpMd = block + "@" + self._getExtraPath('output.xmd')
                iterator = md.SetMdIterator(xmpMd, sortByLabel=md.MDL_ITEM_ID,
                                            updateItemCallback=self._updateParticle,
                                            skipDisabled=True)
                # itemDataIterator is not neccesary because, the class SetMdIterator
                # contain all the information about the metadata
                clsSet.classifyItems(updateItemCallback=iterator.updateItem,
                                     updateClassCallback=self._updateClass)