# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *              Tomas Majtner (tmajtner@cnb.csic.es)  -- streaming version
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
import convert
import xmipp
from pyworkflow.em.data import SetOfCTF
from pyworkflow.object import Set, Float
import pyworkflow.protocol.constants as cons
import pyworkflow.protocol.params as params
import pyworkflow.em as em
from pyworkflow.em.metadata import Row, MetaData


class XmippProtCTFDiscrepancy(em.ProtCTFMicrographs):
    """
    Protocol to estimate the agreement between different estimation of the CTF
    for the same set of micrographs. The algorithm assumes that two CTF are
    consistent if the phase (wave aberration function) of the two CTFs are
    closer than 90 degrees. The reported resolution is the resolution at
    which the two CTF phases differ in 90 degrees.
    """
    _label = 'ctf consensus'

    def __init__(self, **args):
        em.ProtCTFMicrographs.__init__(self, **args)
        self._freqResol = {}
        self.stepsExecutionMode = params.STEPS_SERIAL

    def _defineParams(self, form):
        form.addSection(label='Input')
        # Read N ctfs estimations
        form.addParam('inputCTF1', params.PointerParam,
                      pointerClass='SetOfCTF',
                      label="Reference CTF",
                      help='Reference CTF in comparison')
        form.addParam('inputCTF2', params.PointerParam,
                      pointerClass='SetOfCTF',
                      label="target CTF",
                      help='CTF to be compared with reference CTF')
        form.addParallelSection(threads=0, mpi=0)

# --------------------------- INSERT steps functions -------------------------
    def _insertAllSteps(self):
        self.finished = False
        self.insertedDict = {}
        self.processedDict = []
        ctfSteps = self._checkNewInput()
        self._insertFunctionStep('createOutputStep',
                                 prerequisites=ctfSteps, wait=True)

    def createOutputStep(self):
        pass

    def _getFirstJoinStepName(self):
        # This function will be used for streaming, to check which is
        # the first function that need to wait for all ctfs
        # to have completed, this can be overriden in subclasses
        # (e.g., in Xmipp 'sortPSDStep')
        return 'createOutputStep'

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None

    def _insertNewCtfsSteps(self, insertedDict, SetOfCtf1, SetOfCtf2):
        deps = []
        stepId = self._insertFunctionStep("computeCTFDiscrepancyStep",
                                          SetOfCtf1, SetOfCtf2,
                                          prerequisites=[])
        deps.append(stepId)
        if len(SetOfCtf1) < len(SetOfCtf2):
            for ctf in SetOfCtf1:
                if ctf.getObjId() not in insertedDict:
                    insertedDict[ctf.getObjId()] = stepId
        else:
            for ctf in SetOfCtf2:
                if ctf.getObjId() not in insertedDict:
                    insertedDict[ctf.getObjId()] = stepId
        return deps

    def _stepsCheck(self):
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        # Check if there are new ctf to process from the input set
        ctfsFile1 = self.inputCTF1.get().getFileName()
        ctfsFile2 = self.inputCTF2.get().getFileName()
        ctfsSet1 = SetOfCTF(filename=ctfsFile1)
        ctfsSet2 = SetOfCTF(filename=ctfsFile2)
        ctfsSet1.loadAllProperties()
        ctfsSet2.loadAllProperties()
        self.SetOfCtf1 = [m.clone() for m in ctfsSet1]
        self.SetOfCtf2 = [m.clone() for m in ctfsSet2]
        self.streamClosed = ctfsSet1.isStreamClosed() \
                            and ctfsSet2.isStreamClosed()
        ctfsSet1.close()
        ctfsSet2.close()
        if len(self.SetOfCtf1) < len(self.SetOfCtf2):
            ctfsSet = self.SetOfCtf1
        else:
            ctfsSet = self.SetOfCtf2
        newCtfs = any(m.getObjId() not in self.insertedDict
                      for m in ctfsSet)
        outputStep = self._getFirstJoinStep()
        if newCtfs:
            fDeps = self._insertNewCtfsSteps(self.insertedDict,
                                             ctfsSet1, ctfsSet2)
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()


    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return
        # Load previously done items (from text file)
        doneList = self._readDoneList()
        # Check for newly done items
        if len(self.SetOfCtf1) < len(self.SetOfCtf2):
            SetOfCtf = self.SetOfCtf1
        else:
            SetOfCtf = self.SetOfCtf2
        newDone = [m.clone() for m in SetOfCtf
                   if int(m.getObjId()) not in doneList]
        self.finished = self.streamClosed and \
                        (len(doneList) == len(SetOfCtf))
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


    def _loadOutputSet(self, SetClass, baseName):
        setFile = self._getPath(baseName)
        if os.path.exists(setFile):
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        outputSet.copyInfo(self.SetOfCtf1)
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
            ctfs = self._createSetOfCTF()
            inputCTF = self.inputCTF1.get()
            ctfs.copyInfo(inputCTF)
            ctfs.setMicrographs(inputCTF.getMicrographs())
            for ctf in inputCTF:
                if ctf.getObjId() > len(self.inputCTF2.get()):
                    break
                ctfAux = ctf.clone()
                ctfId = ctf.getObjId()
                resolution = self._freqResol[ctfId]
                ctfAux.setResolution(resolution)
                ctfAux._discrepancy_astigmatism = Float(ctf.getDefocusU() -
                                                        ctf.getDefocusV())
                ctfs.append(ctfAux)
            self._defineOutputs(outputCTF=ctfs)
            self._defineSourceRelation(inputCTF, ctfs)

    def _ctfToMd(self, ctf, ctfMd):
        """ Write the proper metadata for Xmipp from a given CTF """
        ctfMd.clear()
        ctfRow = Row()
        convert.ctfModelToRow(ctf, ctfRow)
        convert.micrographToRow(ctf.getMicrograph(), ctfRow,
                                alignType=convert.ALIGN_NONE)
        ctfRow.addToMd(ctfMd)

    def computeCTFDiscrepancyStep(self, SetOfCtf1, SetOfCtf2):
        self._computeCTFDiscrepancy(SetOfCtf1, SetOfCtf2)

    def _computeCTFDiscrepancy(self, method1, method2):
        # TODO must be same micrographs
        # move to a single step, each step takes 5 sec while the function
        # takes 0.03 sec
        # convert to md
        md1 = MetaData()
        md2 = MetaData()

        for ctf1 in method1:  # reference CTF
            ctfId = ctf1.getObjId()
            if ctfId in self.processedDict or ctfId > len(method2):
                continue
            self.processedDict.append(ctfId)
            self._ctfToMd(ctf1, md1)
            ctf2 = method2[ctfId]  # .target
            self._ctfToMd(ctf2, md2)
            self._freqResol[ctfId] = xmipp.errorMaxFreqCTFs2D(md1, md2)

        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN
        outSet = self._loadOutputSet(SetOfCTF, 'ctfs.sqlite')
        self._updateOutputSet('outputCtfs', outSet, streamMode)

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
        """ Write to a text file the items that have been done. """
        with open(self._getAllDone(), 'a') as f:
            for part in partList:
                f.write('%d\n' % part.getObjId())

    def _citations(self):
        return ['Marabini2014a']

    def _summary(self):
        message = []
        for i, ctfs in enumerate([self.inputCTF1, self.inputCTF2]):
            protocol = self.getMapper().getParent(ctfs.get())
            message.append("Method %d: %s" % (i+1, protocol.getClassLabel()))
        return message

    def _validate(self):
        """ The function of this hook is to add some validation before the
        protocol is launched to be executed. It should return a list of
        errors. If the list is empty the protocol can be executed.
        """
        # same micrographs in both CTF??
        errors = []
        # Add some errors if input is not valid
        return errors