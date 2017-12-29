# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Amaya Jimenez    (ajimenez@cnb.csic.es)
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


from os.path import getmtime, exists
from datetime import datetime

import pyworkflow.em as em
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils

from pyworkflow.object import Set
from pyworkflow.protocol.constants import (STATUS_NEW)


class XmippProtCTFSelection(em.ProtCTFMicrographs):
    """
    Protocol to make a selection of meaningful CTFs in basis of the defocus
    values, the astigmatism, and the resolution.
    """
    _label = 'ctf selection'

    def __init__(self, **args):
        em.ProtCTFMicrographs.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        # Read a ctf estimation
        form.addParam('inputCTF', params.PointerParam,
                      pointerClass='SetOfCTF',
                      label="Input CTF",
                      important=True,
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

    # --------------------------- INSERT steps functions ------------------
    def _insertAllSteps(self):

        fDeps = []
        self.insertedDict = {}

        fDeps = self._insertNewSelectionSteps(self.insertedDict,
                                              self.inputCTF.get())
        # For the streaming mode, the steps function have a 'wait' flag
        # that can be turned on/off. For example, here we insert the
        # createOutputStep but it wait=True, which means that can not be
        # executed until it is set to False
        # (when the input micrographs stream is closed)
        self._insertFunctionStep('createOutputStep', prerequisites=fDeps,
                                 wait=True)

    def _insertNewSelectionSteps(self, insertedDict, inputCtfs):
        """ Insert steps to process new ctfs (from streaming)
        Params:
            insertedDict: contains already processed ctfs
            inputCtfs: input ctfs set to be check
        """
        deps = []
        # For each ctf insert the step to process it
        for ctf in inputCtfs:
            ctfId = ctf.getObjId()
            if ctfId not in insertedDict:
                stepId = self._insertCtfSelectionStep(ctfId)
                deps.append(stepId)
                insertedDict[ctfId] = stepId
        return deps

    def _insertCtfSelectionStep(self, ctfId):
        """ Insert the processMovieStep for a given movie. """
        # Note1: At this point is safe to pass the movie, since this
        # is not executed in parallel, here we get the params
        # to pass to the actual step that is gone to be executed later on
        # Note2: We are serializing the Movie as a dict that can be passed
        # as parameter for a functionStep
        stepId = self._insertFunctionStep('_selectCTF', ctfId,
                                          prerequisites=[])
        return stepId

    # --------------------------- STEPS functions ------------------------
    def _stepsCheck(self):

        # check if there are new ctfs and process them
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        """ Check if there are new ctf to be processed and add the necessary
        steps."""
        ctfFile = self.inputCTF.get().getFileName()

        now = datetime.now()
        self.lastCheck = getattr(self, 'lastCheck', now)
        mTime = datetime.fromtimestamp(getmtime(ctfFile))
        self.debug('Last check: %s, modification: %s'
                   % (pwutils.prettyTime(self.lastCheck),
                      pwutils.prettyTime(mTime)))

        # Open input ctfs.sqlite and close it as soon as possible
        self._loadInputList()
        # If the input ctfs.sqlite have not changed since our last check,
        # it does not make sense to check for new input data
        if self.lastCheck > mTime and hasattr(self, 'listOfCtf'):
            return None

        self.lastCheck = now
        newCtf = any(ctf.getObjId() not in
                     self.insertedDict for ctf in self.listOfCtf)
        outputStep = self._getFirstJoinStep()

        if newCtf:
            fDeps = self._insertNewSelectionSteps(self.insertedDict,
                                                  self.listOfCtf)
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
        self.finished = (self.isStreamClosed == Set.STREAM_CLOSED and
                         allDone == len(self.listOfCtf))
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

                ctfSet.append(ctf)
                micSet.append(mic)
                self._writeDoneListAccepted(ctfId)

            inputCtfSet.close()

        if newDoneDiscarded:
            inputCtfSet = self._loadInputCtfSet()
            for ctfId in newDoneDiscarded:
                ctf = inputCtfSet[ctfId].clone()
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

        if (exists(self._getPath('ctfs.sqlite'))):
            self._updateOutputSet('outputCTF', ctfSet, streamMode)
            self._updateOutputSet('outputMicrographs', micSet, streamMode)
        # AJ new subsets with discarded ctfs
        if (exists(self._getPath('ctfsDiscarded.sqlite'))):
            self._updateOutputSet('outputCTFDiscarded',
                                  ctfSetDiscarded, streamMode)
            self._updateOutputSet('outputMicrographsDiscarded',
                                  micSetDiscarded, streamMode)

        if (exists(self._getPath('ctfs.sqlite'))):
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
        if (exists(self._getPath('ctfsDiscarded.sqlite'))):
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

        if (exists(self._getPath('ctfs.sqlite'))):
            ctfSet.close()
            micSet.close()
        # AJ new subsets with discarded ctfs
        if (exists(self._getPath('ctfsDiscarded.sqlite'))):
            micSetDiscarded.close()
            ctfSetDiscarded.close()

    def createOutputStep(self):
        # Do nothing now, the output should be ready.
        pass

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

    # --------------------------- INFO functions --------------------------
    def _summary(self):
        message = []
        return message

    def _methods(self):
        # nothing here
        pass

    def _validate(self):
        """ The function of this hook is to add some validation before the
        protocol is launched to be executed. It should return a list of
        errors. If the list is empty the protocol can be executed.
        """
        message = []
        return message

    # --------------------------- UTILS functions ---------------------

    def usingXmipp(self, ctf):
        return ctf.hasAttribute('_xmipp_ctfCritfirstZeroRatio')

    def _getCtfResol(self, ctf):
        resolution = ctf.getResolution()
        if resolution is not None:
            return resolution
        else:
            return 0
        # if ctf.hasAttribute('_ctffind4_ctfResolution'):
        #     return ctf._ctffind4_ctfResolution.get(), False
        # elif ctf.hasAttribute('_gctf_ctfResolution'):
        #     return ctf._gctf_ctfResolution.get(), False
        # elif ctf.hasAttribute('_xmipp_ctfCritMaxFreq'):
        #     return ctf._xmipp_ctfCritMaxFreq.get(), True
        # else:
        #     # if 0, the protocol does not select by ctf resolution.
        #     return 0

    def _readDoneListDiscarded(self):
        """ Read from a text file the id's of the items
        that have been done. """
        DiscardedFile = self._getDiscardedDone()
        DiscardedList = []
        # Check what items have been previously done
        if exists(DiscardedFile):
            with open(DiscardedFile) as f:
                DiscardedList += [int(line.strip()) for line in f]

        return DiscardedList

    def _readDoneListAccepted(self):
        """ Read from a text file the id's of the items
        that have been done. """
        AcceptedFile = self._getAcceptedDone()
        AcceptedList = []
        # Check what items have been previously done
        if exists(AcceptedFile):
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

    def _loadInputList(self):
        """ Load the input set of ctfs and create a list. """
        ctfSet = self._loadInputCtfSet()
        self.isStreamClosed = ctfSet.getStreamState()
        self.listOfCtf = [m.clone() for m in ctfSet]
        ctfSet.close()
        self.debug("Closed db.")

    def _loadOutputSet(self, SetClass, baseName):
        """
        Load the output set if it exists or create a new one.
        """
        setFile = self._getPath(baseName)

        if exists(setFile):
            outputSet = SetClass(filename=setFile)
            if(outputSet.__len__() is 0):
                pwutils.path.cleanPath(setFile)

        if exists(setFile):
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

    def _readtCtfId(self, accepted):
        if accepted:
            fn = self._getCtfSelecFileAccepted()
        else:
            fn = self._getCtfSelecFileDiscarded()
        ctfList = []
        # Check what items have been previously done
        if exists(fn):
            with open(fn) as f:
                ctfList += [int(line.strip().split()[0]) for line in f]
        return ctfList

    def _getEnable(self, ctfId):
        fn = self._getCtfSelecFileAccepted()
        ctfList = []
        # Check what items have been previously done
        if exists(fn):
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
