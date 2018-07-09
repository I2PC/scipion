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
        return self._insertNewMics(inputMics,
                                   lambda mic: mic.getMicName(),
                                   self._insertExtractMicrographStep,
                                   self._insertExtractMicrographListStep,
                                   *self._getExtractArgs())

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
    
    # -------------------------- STEPS functions ------------------------------


    def extractMicrographStep(self, micKey, *args):
        """ Step function that will be common for all extraction protocols.
        It will take an id and will grab the micrograph from a micDict map.
        The micrograph will be passed as input to the _extractMicrograph
        function.
        """
        # Retrieve the corresponding micrograph with this key and the
        # associated list of coordinates
        mic = self.micDict[micKey]

        micDoneFn = self._getMicDone(mic)
        micFn = mic.getFileName()

        if self.isContinued() and os.path.exists(micDoneFn):
            self.info("Skipping micrograph: %s, seems to be done" % micFn)
            return

        coordList = self.coordDict[mic.getObjId()]
        self._convertCoordinates(mic, coordList)

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

    # ---------- Methods to extract many micrographs at once ------------------

    def _insertExtractMicrographListStep(self, micList, prerequisites, *args):
        """ Basic method to insert a picking step for a given micrograph. """
        return self._insertFunctionStep('extractMicrographListStep',
                                        [mic.getMicName() for mic in micList],
                                        *args, prerequisites=prerequisites)

    def extractMicrographListStep(self, micKeyList, *args):
        micList = []

        for micName in micKeyList:
            mic = self.micDict[micName]
            micDoneFn = self._getMicDone(mic)
            micFn = mic.getFileName()
            if self.isContinued() and os.path.exists(micDoneFn):
                self.info("Skipping micrograph: %s, seems to be done" % micFn)

            else:
                # Clean old finished files
                pwutils.cleanPath(micDoneFn)
                self.info("Extracting micrograph: %s " % micFn)
                micList.append(mic)

        self._extractMicrographList(micList, *args)

        for mic in micList:
            # Mark this mic as finished
            open(self._getMicDone(mic), 'w').close()

    def _extractMicrographList(self, micList, *args):
        """ Extract more than one micrograph at once.
        Here the default implementation is to iterate through the list and
        call the single extract, but it could be re-implemented on each
        subclass to provide a more efficient implementation.
        """
        for mic in micList:
            self._extractMicrograph(mic, *args)

    # --------------------------- UTILS functions -----------------------------
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

        # Store this value to be used when inserting new steps and batch mode
        self.streamClosed = self._isStreamClosed()

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

            # Maybe it would be good idea to take a snap to avoid
            # so much IO if this protocol does not have much to do now
            if allDone == len(self.micDict):
                self._streamingSleepOnWait()

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

    def readPartsFromMics(self, micDoneList, outputParts):
        """ This method should be implemented in subclasses to read
        the particles from a given list of micrographs.
        """
        pass

    def _updateOutputPartSet(self, micList, streamMode):
        outputName = 'outputParticles'
        outputParts = getattr(self, outputName, None)
        firstTime = True

        if outputParts is None:
            inputMics = self.getInputMicrographs()
            outputParts = self._createSetOfParticles()
            outputParts.copyInfo(inputMics)
            outputParts.setCoordinates(self.getCoords())

            if self.getAttributeValue('doFlip', False):
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
