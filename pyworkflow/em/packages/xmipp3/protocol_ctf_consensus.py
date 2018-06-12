# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *              Tomas Majtner (tmajtner@cnb.csic.es)  -- streaming version
# *              Amaya Jimenez (ajimenez@cnb.csic.es)
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
from datetime import datetime
from pyworkflow.em.data import SetOfCTF
from pyworkflow.object import Set, Float
import pyworkflow.protocol.params as params
import pyworkflow.em as em
import pyworkflow.utils as pwutils
from pyworkflow.em.metadata import Row, MetaData
from pyworkflow.protocol.constants import (STATUS_NEW)


class XmippProtCTFConsensus(em.ProtCTFMicrographs):
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
        form.addParam('inputCTF', params.PointerParam, pointerClass='SetOfCTF',
                      label="Input CTF", important=True,
                      help='Select the estimated CTF to evaluate')
        form.addParam('useDefocus', params.BooleanParam, default=True,
                      label='Use Defocus for selection',
                      help='Use this button to decide if carry out the '
                           'selection taking into account or not the defocus '
                           'values.')

        line = form.addLine('Defocus (A)', condition="useDefocus",
                            help='Minimum and maximum values for defocus in '
                                 'Angstroms')
        line.addParam('minDefocus', params.FloatParam, default=0, label='Min')
        line.addParam('maxDefocus', params.FloatParam,
                      default=40000, label='Max')

        form.addParam('useAstigmatism', params.BooleanParam, default=True,
                      label='Use Astigmatism for selection',
                      help='Use this button to decide if carry out the '
                           'selection taking into account or not the '
                           'astigmatism value.')
        form.addParam('astigmatism', params.FloatParam, default=1000,
                      label='Astigmatism (A)', condition="useAstigmatism",
                      help='Maximum value allowed for astigmatism in '
                           'Angstroms. If the evaluated CTF does'
                           ' not fulfill '
                           'this requirement, it will be discarded.')
        form.addParam('useResolution', params.BooleanParam, default=True,
                      label='Use Resolution for selection',
                      help='Use this button to decide if carry out the '
                           'selection taking into account or not the '
                           'resolution value.')
        form.addParam('resolution', params.FloatParam, default=7,
                      label='Resolution (A)',
                      condition="useResolution",
                      help='Minimum value for resolution in Angstroms. '
                           'If the evaluated CTF does not fulfill this '
                           'requirement, it will be discarded.')

        form.addParam('useCritXmipp', params.BooleanParam, default=False,
                      label='Use Xmipp Crit criterion for selection',
                      help='Use this button to decide if carrying out the '
                           'selection taking into account the Xmipp '
                           'Crit parameters. \n'
                           'Only available when Xmipp CTF'
                           ' estimation was used.')
        form.addParam('critFirstZero', params.FloatParam, default=5,
                      condition="useCritXmipp", label='CritFirstZero',
                      help='Minimun value of CritFirstZero')
        line = form.addLine('CritFirstZeroRatio', condition="useCritXmipp",
                            help='Minimum and maximum values for '
                                 'CritFirstZeroRatio')
        line.addParam('minCritFirstZeroRatio', params.FloatParam, default=0.9,
                      label='Min')
        line.addParam('maxCritFirstZeroRatio', params.FloatParam, default=1.1,
                      label='Max')
        form.addParam('critCorr', params.FloatParam, default=0,
                      condition="useCritXmipp", label='CritCorr',
                      help='Minimum value of CritCorr')
        form.addParam('critCtfMargin', params.FloatParam, default=0,
                      condition="useCritXmipp", label='CritCtfMargin',
                      help='Minimum value of CritCtfMargin')
        line = form.addLine('CritNonAstigmaticValidity',
                            condition="useCritXmipp",
                            help='Minimum and maximum values for '
                                 'CritNonAstigmaticValidity')
        line.addParam('minCritNonAstigmaticValidity', params.FloatParam,
                      default=0.3, label='Min')
        line.addParam('maxCritNonAstigmaticValidity', params.FloatParam,
                      default=25, label='Max')
        form.addParam('calculateConsensus', params.BooleanParam, default=False,
                      label='Calculate Consensus Resolution',
                      help='Option for calculating consensus resolution.')
        form.addParam('inputCTF2', params.PointerParam,
                      pointerClass='SetOfCTF', condition="calculateConsensus",
                      label="target CTF",
                      help='CTF to be compared with reference CTF')
        form.addParam('minConsResol', params.FloatParam,
                      condition="calculateConsensus", default=5.0,
                      label='Minimum consensus resolution in angstroms.')
        form.addParallelSection(threads=0, mpi=0)

# --------------------------- INSERT steps functions -------------------------
    def _insertAllSteps(self):
        self.finished = False
        self.insertedDict = {}
        self.processedDict = []
        self.outputDict = []
        self.allCtf1 = []
        self.allCtf2 = []
        if self.calculateConsensus:
            ctfSteps = self._checkNewInput()
        else:
            ctfSteps = self._insertNewSelectionSteps(self.insertedDict,
                                                     self.inputCTF.get())
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

    def _insertNewCtfsSteps(self, SetOfCtf1, SetOfCtf2, insDict):
        deps = []
        discrepId = self._insertFunctionStep("_computeCTFDiscrepancy",
                                             SetOfCtf1, SetOfCtf2,
                                             prerequisites=[])
        deps.append(discrepId)
        if (len(SetOfCtf1) > len(SetOfCtf2)):
            SetOfCtf = SetOfCtf2
        else:
            SetOfCtf = SetOfCtf1
        for ctf in SetOfCtf:
            ctfId = ctf.getObjId()
            if ctfId not in insDict:
                stepId = self._insertFunctionStep('_selectCTF', ctfId,
                                                  prerequisites=[discrepId])
                deps.append(stepId)
                insDict[ctfId] = stepId

        return deps

    def _insertNewSelectionSteps(self, insertedDict, inputCtfs):
        deps = []
        # For each ctf insert the step to process it
        for ctf in inputCtfs:
            ctfId = ctf.getObjId()
            if ctfId not in insertedDict:
                stepId = self._insertFunctionStep('_selectCTF', ctfId,
                                                  prerequisites=[])
                deps.append(stepId)
                insertedDict[ctfId] = stepId
        return deps

    def _stepsCheck(self):
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        if self.calculateConsensus:
            # Check if there are new ctf to process from the input set
            ctfsFile1 = self.inputCTF.get().getFileName()
            ctfsFile2 = self.inputCTF2.get().getFileName()
            self.lastCheck = getattr(self, 'lastCheck', datetime.now())
            mTime = max(datetime.fromtimestamp(os.path.getmtime(ctfsFile1)),
                        datetime.fromtimestamp(os.path.getmtime(ctfsFile2)))
            # If the input movies.sqlite have not changed since our last check,
            # it does not make sense to check for new input data
            if self.lastCheck > mTime and hasattr(self, 'SetOfCtf1'):
                return None
            ctfsSet1 = SetOfCTF(filename=ctfsFile1)
            ctfsSet2 = SetOfCTF(filename=ctfsFile2)
            ctfsSet1.loadAllProperties()
            ctfsSet2.loadAllProperties()
            if len(self.allCtf1) > 0:
                newCtf1 = [ctf.clone() for ctf in
                           ctfsSet1.iterItems(orderBy='creation',
                                              where='creation>"' + str(
                                                  self.checkCtf1) + '"')]
            else:
                newCtf1 = [ctf.clone() for ctf in ctfsSet1]
            self.allCtf1 = self.allCtf1 + newCtf1
            if len(newCtf1) > 0:
                for ctf in ctfsSet1.iterItems(orderBy='creation',
                                              direction='DESC'):
                    self.checkCtf1 = ctf.getObjCreation()
                    break
            if len(self.allCtf2) > 0:
                newCtf2 = [ctf.clone() for ctf in
                           ctfsSet2.iterItems(orderBy='creation',
                                              where='creation>"' + str(
                                                  self.checkCtf2) + '"')]
            else:
                newCtf2 = [ctf.clone() for ctf in ctfsSet2]
            self.allCtf2 = self.allCtf2 + newCtf2
            if len(newCtf2) > 0:
                for ctf in ctfsSet2.iterItems(orderBy='creation',
                                              direction='DESC'):
                    self.checkCtf2 = ctf.getObjCreation()
                    break
            self.lastCheck = datetime.now()
            self.isStreamClosed = ctfsSet1.isStreamClosed() and \
                                  ctfsSet2.isStreamClosed()
            ctfsSet1.close()
            ctfsSet2.close()
            outputStep = self._getFirstJoinStep()
            if len(set(self.allCtf1)) > len(set(self.processedDict)) and \
               len(set(self.allCtf2)) > len(set(self.processedDict)):
                fDeps = self._insertNewCtfsSteps(ctfsSet1, ctfsSet2,
                                                 self.insertedDict)
                if outputStep is not None:
                    outputStep.addPrerequisites(*fDeps)
                self.updateSteps()
        else:
            ctfFile = self.inputCTF.get().getFileName()
            now = datetime.now()
            self.lastCheck = getattr(self, 'lastCheck', now)
            mTime = datetime.fromtimestamp(os.path.getmtime(ctfFile))
            self.debug('Last check: %s, modification: %s'
                       % (pwutils.prettyTime(self.lastCheck),
                          pwutils.prettyTime(mTime)))

            # Open input ctfs.sqlite and close it as soon as possible
            ctfSet = self._loadInputCtfSet()
            self.isStreamClosed = ctfSet.isStreamClosed()
            self.allCtf1 = [m.clone() for m in ctfSet]
            ctfSet.close()

            # If the input ctfs.sqlite have not changed since our last check,
            # it does not make sense to check for new input data
            if self.lastCheck > mTime and hasattr(self, 'allCtf1'):
                return None

            self.lastCheck = now
            newCtf = any(ctf.getObjId() not in
                         self.insertedDict for ctf in self.allCtf1)
            outputStep = self._getFirstJoinStep()

            if newCtf:
                fDeps = self._insertNewSelectionSteps(self.insertedDict,
                                                      self.allCtf1)
                if outputStep is not None:
                    outputStep.addPrerequisites(*fDeps)
                self.updateSteps()


    def _checkNewOutput(self):
        """ Check for already selected CTF and update the output set. """

        # Load previously done items (from text file)
        doneListDiscarded = self._readDoneListDiscarded()
        doneListAccepted = self._readDoneListAccepted()

        # Check for newly done items
        ctfListIdAccepted = self._readtCtfId(True)
        ctfListIdDiscarded = self._readtCtfId(False)

        newDoneAccepted = [ctfId for ctfId in ctfListIdAccepted
                           if ctfId not in doneListAccepted]
        newDoneDiscarded = [ctfId for ctfId in ctfListIdDiscarded
                            if ctfId not in doneListDiscarded]
        firstTimeAccepted = len(doneListAccepted) == 0
        firstTimeDiscarded = len(doneListDiscarded) == 0
        allDone = len(doneListAccepted) + len(doneListDiscarded) +\
                  len(newDoneAccepted) + len(newDoneDiscarded)

        # We have finished when there is not more input ctf (stream closed)
        # and the number of processed ctf is equal to the number of inputs
        self.finished = (self.isStreamClosed and allDone == len(self.allCtf1))

        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        # reading the outputs
        if (len(doneListAccepted) > 0 or len(newDoneAccepted) > 0):
            ctfSet = self._loadOutputSet(em.SetOfCTF, 'ctfs.sqlite')
            micSet = self._loadOutputSet(em.SetOfMicrographs,
                                         'micrographs.sqlite')

        # AJ new subsets with discarded ctfs
        if (len(doneListDiscarded) > 0 or len(newDoneDiscarded) > 0):
            ctfSetDiscarded = \
                self._loadOutputSet(em.SetOfCTF, 'ctfsDiscarded.sqlite')
            micSetDiscarded = \
                self._loadOutputSet(em.SetOfMicrographs,
                                    'micrographsDiscarded.sqlite')

        if newDoneAccepted:
            inputCtfSet = self._loadInputCtfSet()
            for ctfId in newDoneAccepted:
                ctf = inputCtfSet[ctfId].clone()
                mic = ctf.getMicrograph().clone()

                ctf.setEnabled(self._getEnable(ctfId))
                mic.setEnabled(self._getEnable(ctfId))

                if self.calculateConsensus:
                    ctf.setConsensusResolution(self._freqResol[ctfId])
                    ctf._discrepancy_astigmatism = Float(ctf.getDefocusU() -
                                                        ctf.getDefocusV())

                ctfSet.append(ctf)
                micSet.append(mic)
                self._writeDoneListAccepted(ctfId)

            inputCtfSet.close()

        if newDoneDiscarded:
            inputCtfSet = self._loadInputCtfSet()
            for ctfId in newDoneDiscarded:
                ctf = inputCtfSet[ctfId].clone()
                if self.calculateConsensus:
                    ctf.setConsensusResolution(self._freqResol[ctfId])
                    ctf._discrepancy_astigmatism = Float(ctf.getDefocusU() -
                                                        ctf.getDefocusV())
                mic = ctf.getMicrograph().clone()
                micSetDiscarded.append(mic)
                ctfSetDiscarded.append(ctf)
                self._writeDoneListDiscarded(ctfId)

            inputCtfSet.close()

        if not self.finished and not newDoneDiscarded and not newDoneAccepted:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        if (os.path.exists(self._getPath('ctfs.sqlite'))):
            self._updateOutputSet('outputCTF', ctfSet, streamMode)
            self._updateOutputSet('outputMicrographs', micSet, streamMode)
        # AJ new subsets with discarded ctfs
        if (os.path.exists(self._getPath('ctfsDiscarded.sqlite'))):
            self._updateOutputSet('outputCTFDiscarded',
                                  ctfSetDiscarded, streamMode)
            self._updateOutputSet('outputMicrographsDiscarded',
                                  micSetDiscarded, streamMode)

        if (os.path.exists(self._getPath('ctfs.sqlite'))):
            if firstTimeAccepted:
                # define relation just once
                self._defineSourceRelation(
                    self.inputCTF.get().getMicrographs(), micSet)
                self._defineSourceRelation(ctfSet, micSet)
                self._defineSourceRelation(self.inputCTF, ctfSet)
                self._defineCtfRelation(micSet, ctfSet)
            else:
                ctfSet.close()
                micSet.close()

        # AJ new subsets with discarded ctfs
        if (os.path.exists(self._getPath('ctfsDiscarded.sqlite'))):
            if firstTimeDiscarded:
                self._defineSourceRelation(
                    self.inputCTF.get().getMicrographs(), micSetDiscarded)
                self._defineSourceRelation(ctfSetDiscarded, micSetDiscarded)
                self._defineSourceRelation(self.inputCTF, ctfSetDiscarded)
                self._defineCtfRelation(micSetDiscarded, ctfSetDiscarded)
            else:
                micSetDiscarded.close()
                ctfSetDiscarded.close()

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)

        if (os.path.exists(self._getPath('ctfs.sqlite'))):
            ctfSet.close()
            micSet.close()
        # AJ new subsets with discarded ctfs
        if (os.path.exists(self._getPath('ctfsDiscarded.sqlite'))):
            micSetDiscarded.close()
            ctfSetDiscarded.close()


    def _loadOutputSet(self, SetClass, baseName):
        """
        Load the output set if it exists or create a new one.
        """
        setFile = self._getPath(baseName)

        if os.path.exists(setFile):
            outputSet = SetClass(filename=setFile)
            if (outputSet.__len__() is 0):
                pwutils.path.cleanPath(setFile)

        if os.path.exists(setFile):
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        micSet = self.inputCTF.get().getMicrographs()

        if isinstance(outputSet, em.SetOfMicrographs):
            outputSet.copyInfo(micSet)
        elif isinstance(outputSet, em.SetOfCTF):
            outputSet.setMicrographs(micSet)
        return outputSet

    def _ctfToMd(self, ctf, ctfMd):
        """ Write the proper metadata for Xmipp from a given CTF """
        ctfMd.clear()
        ctfRow = Row()
        convert.ctfModelToRow(ctf, ctfRow)
        convert.micrographToRow(ctf.getMicrograph(), ctfRow,
                                alignType=convert.ALIGN_NONE)
        ctfRow.addToMd(ctfMd)

    def _computeCTFDiscrepancy(self, method1, method2):
        # TODO must be same micrographs
        # move to a single step, each step takes 5 sec while the function
        # takes 0.03 sec
        # convert to md
        md1 = MetaData()
        md2 = MetaData()

        for ctf1 in method1:  # reference CTF
            ctfId = ctf1.getObjId()
            if ctfId in self.processedDict:
                continue
            for ctf2 in method2:
                ctfId2 = ctf2.getObjId()
                if ctfId2 != ctfId:
                    continue
                self.processedDict.append(ctfId)
                self._ctfToMd(ctf1, md1)
                self._ctfToMd(ctf2, md2)
                self._freqResol[ctfId] = xmipp.errorMaxFreqCTFs2D(md1, md2)

    def _selectCTF(self, ctfId):
        # Depending on the flags selected by the user, we set the values of
        # the params to compare with

        minDef, maxDef = self._getDefociValues()
        maxAstig = self._getMaxAstisgmatism()
        minResol = self._getMinResol()

        # TODO: Change this way to get the ctf.
        ctf = self.inputCTF.get()[ctfId]

        defocusU = ctf.getDefocusU()
        defocusV = ctf.getDefocusV()
        astigm = defocusU - defocusV
        resol = self._getCtfResol(ctf)
        usingXmipp = self.usingXmipp(ctf)

        firstCondition = ((defocusU < minDef) or (defocusU > maxDef) or
                          (defocusV < minDef) or (defocusV > maxDef) or
                          (astigm > maxAstig) or (resol > minResol))

        if self.calculateConsensus:
            firstCondition = firstCondition or (
            self.minConsResol < self._freqResol[ctfId])

        secondCondition = False
        if (self.useCritXmipp and usingXmipp):
            firstZero = self._getCritFirstZero()
            minFirstZero, maxFirstZero = self._getCritFirstZeroRatio()
            corr = self._getCritCorr()
            ctfMargin = self._getCritCtfMargin()
            minNonAstigmatic, maxNonAstigmatic = \
                self._getCritNonAstigmaticValidity()

            secondCondition = (
                (ctf._xmipp_ctfCritFirstZero.get() < firstZero) or
                (ctf._xmipp_ctfCritfirstZeroRatio.get() < minFirstZero) or
                (ctf._xmipp_ctfCritfirstZeroRatio.get() > maxFirstZero) or
                (ctf._xmipp_ctfCritCorr13.get() < corr) or
                (ctf._xmipp_ctfCritCtfMargin.get() < ctfMargin) or
                (ctf._xmipp_ctfCritNonAstigmaticValidty.get() <
                 minNonAstigmatic) or
                (ctf._xmipp_ctfCritNonAstigmaticValidty.get() >
                 maxNonAstigmatic))

        """ Write to a text file the items that have been done. """
        if firstCondition or secondCondition:
            fn = self._getCtfSelecFileDiscarded()
            with open(fn, 'a') as f:
                f.write('%d F\n' % ctf.getObjId())
        else:
            if (ctf.isEnabled()):
                fn = self._getCtfSelecFileAccepted()
                with open(fn, 'a') as f:
                    f.write('%d T\n' % ctf.getObjId())
            else:
                fn = self._getCtfSelecFileAccepted()
                with open(fn, 'a') as f:
                    f.write('%d F\n' % ctf.getObjId())

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
        for i, ctfs in enumerate([self.inputCTF, self.inputCTF2]):
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

    def usingXmipp(self, ctf):
        return ctf.hasAttribute('_xmipp_ctfCritfirstZeroRatio')

    def _getCtfResol(self, ctf):
        resolution = ctf.getResolution()
        if resolution is not None:
            return resolution
        else:
            return 0

    def _readDoneListDiscarded(self):
        """ Read from a text file the id's of the items
        that have been done. """
        DiscardedFile = self._getDiscardedDone()
        DiscardedList = []
        # Check what items have been previously done
        if os.path.exists(DiscardedFile):
            with open(DiscardedFile) as f:
                DiscardedList += [int(line.strip()) for line in f]
        return DiscardedList

    def _readDoneListAccepted(self):
        """ Read from a text file the id's of the items
        that have been done. """
        AcceptedFile = self._getAcceptedDone()
        AcceptedList = []
        # Check what items have been previously done
        if os.path.exists(AcceptedFile):
            with open(AcceptedFile) as f:
                AcceptedList += [int(line.strip()) for line in f]
        return AcceptedList

    def _writeDoneListDiscarded(self, ctfId):
        """ Write to a text file the items that have been done. """
        DiscardedFile = self._getDiscardedDone()
        with open(DiscardedFile, 'a') as f:
            f.write('%d\n' % ctfId)

    def _writeDoneListAccepted(self, ctfId):
        """ Write to a text file the items that have been done. """
        AcceptedFile = self._getAcceptedDone()
        with open(AcceptedFile, 'a') as f:
            f.write('%d\n' % ctfId)

    def _getDiscardedDone(self):
        return self._getExtraPath('DONE_discarded.TXT')

    def _getAcceptedDone(self):
        return self._getExtraPath('DONE_accepted.TXT')

    def _getCtfSelecFileAccepted(self):
        return self._getExtraPath('selection-ctf-accepted.txt')

    def _getCtfSelecFileDiscarded(self):
        return self._getExtraPath('selection-ctf-discarded.txt')

    def _readtCtfId(self, accepted):
        if accepted:
            fn = self._getCtfSelecFileAccepted()
        else:
            fn = self._getCtfSelecFileDiscarded()
        ctfList = []
        # Check what items have been previously done
        if os.path.exists(fn):
            with open(fn) as f:
                ctfList += [int(line.strip().split()[0]) for line in f]
        return ctfList

    def _getEnable(self, ctfId):
        fn = self._getCtfSelecFileAccepted()
        # Check what items have been previously done
        if os.path.exists(fn):
            with open(fn) as f:
                for line in f:
                    if ctfId == int(line.strip().split()[0]):
                        if line.strip().split()[1] == 'T':
                            return True
                        else:
                            return False

    def _loadInputCtfSet(self):
        ctfFile = self.inputCTF.get().getFileName()
        self.debug("Loading input db: %s" % ctfFile)
        ctfSet = em.SetOfCTF(filename=ctfFile)
        ctfSet.loadAllProperties()
        return ctfSet

    def _getDefociValues(self):
        if not self.useDefocus:
            return 0, 1000000
        else:
            return self.minDefocus.get(), self.maxDefocus.get()

    def _getMaxAstisgmatism(self):
        if not self.useAstigmatism:
            return 10000000
        else:
            return self.astigmatism.get()

    def _getMinResol(self):
        if not self.useResolution:
            return 1000000
        else:
            return self.resolution.get()

    def _getCritFirstZero(self):
        return self.critFirstZero.get()

    def _getCritFirstZeroRatio(self):
        return self.minCritFirstZeroRatio.get(),\
               self.maxCritFirstZeroRatio.get()

    def _getCritCorr(self):
        return self.critCorr.get()

    def _getCritCtfMargin(self):
        return self.critCtfMargin.get()

    def _getCritNonAstigmaticValidity(self):
        return self.minCritNonAstigmaticValidity.get(), \
               self.maxCritNonAstigmaticValidity.get()
