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


import pyworkflow.protocol.params as params
import pyworkflow.em as em
from datetime import datetime
import os
import pyworkflow.utils as pwutils
from pyworkflow.em.data import SetOfCTF
from pyworkflow.utils.path import copyTree, removeBaseExt, makePath
from pyworkflow.protocol.constants import (STEPS_PARALLEL, LEVEL_ADVANCED,
                                           STATUS_NEW)


class XmippProtCTFSelection(em.ProtCTFMicrographs):
    """
    Protocol to make a selection of meaningful CTFs in basis of the defocus values,
    the astigmatism, and the resolution.
    """
    _label = 'ctf selection'
    
    def __init__(self, **args):
        em.ProtCTFMicrographs.__init__(self, **args)
        #self._freqResol = {}

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        # Read a ctf estimation
        form.addParam('inputCTFs', params.PointerParam, pointerClass='SetOfCTF',
                      important=True,
                      label="Set of CTFs",
                      help='Estimated CTF to evaluate.')
        line = form.addLine('Defocus (A)',
                            help='Minimum and maximum values for defocus in Angstroms')
        line.addParam('minDefocus', params.FloatParam, default=0, label='Min')
        line.addParam('maxDefocus', params.FloatParam, default=40000, label='Max')

        form.addParam('astigmatism', params.FloatParam, default=1000,
                      label='Astigmatism (A)',
                      help='Maximum value allowed for astigmatism in Angstroms. '
                           'If the evaluated CTF does not fulfill this requirement, it will be discarded.')
        form.addParam('resolution', params.FloatParam, default=7,
                  label='Resolution (A)',
                  help='Minimum value for resolution in Angstroms. '
                       'If the evaluated CTF does not fulfill this requirement, it will be discarded.')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """for each ctf insert the steps to compare it
        """
        self.insertedDict = {}
        self.processedMicros = []

        deps = self._insertSteps(self.insertedDict, self.inputCTFs.get())
        fDeps = self._insertFinalSteps(deps)
        waitCondition = True #self._getFirstJoinStepName() == 'createOutputStep'
        self._insertFunctionStep('createOutputStep', prerequisites=fDeps, wait=waitCondition)


    # --------------------------- STEPS functions -------------------------------

    def createOutputStep(self):
        # Do nothing now, the output should be ready.
        print "en createOutputStep"


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


    def readCTF(self):
        print "en readCTF"
        inputCTFs = self.inputCTFs.get()
        outputCTFs = self._createSetOfCTF()

        for ctf in inputCTFs:
            if ctf.getMicrograph().getMicName() not in self.processedMicros:
                self.processedMicros.append(ctf.getMicrograph().getMicName())
                defocusU = ctf.getDefocusU()
                defocusV = ctf.getDefocusV()
                astigm = defocusU - defocusV
                resol = ctf._ctffind4_ctfResolution.get() #PENDIENTEEEEEEEEEEEEEEEEEE
                if defocusU>self.minDefocus and defocusU<self.maxDefocus and \
                    defocusV>self.minDefocus and defocusV<self.maxDefocus and \
                    astigm <self.astigmatism and resol<self.resolution:

                    outputCTFs.append(ctf)

        self._defineOutputs(outputCTF=outputCTFs)
        self._defineTransformRelation(self.inputCTFs.get(), outputCTFs)


    def _insertSteps(self, insertedDict, inputCTFs):
        estimDeps = []
        for ctf in inputCTFs:
            if ctf.getMicrograph().getMicName() not in insertedDict:
                stepId = self._insertFunctionStep('readCTF')
                estimDeps.append(stepId)
                insertedDict[ctf.getMicrograph().getMicName()] = stepId
                break
        return estimDeps


    def _checkNewCTFs(self, ctfSet, outputStep):
        """ Check for new CTF and update the output set. """
        print "en _checkNewCTFs"
        newCTFs = []
        if ctfSet is not None:
            #if ctfSet.get().getFileName() not in self.insertedDict:
            #    newCTFs.append(ctfSet)
            for ctf in ctfSet:
                if ctf.getMicrograph().getMicName() not in self.insertedDict:
                    newCTFs.append(ctfSet)
                    break

        if newCTFs:
            fDeps = self._insertSteps(self.insertedDict, ctfSet)
            self._storeSteps()
            self._numberOfSteps.set(len(self._steps))
            self._store(self._numberOfSteps)
            if outputStep:
                outputStep.addPrerequisites(*fDeps)

        return newCTFs
            #outputCTFs = self._createSetOfCTF()
            #
            #for ctf in ctfSet:
            #    defocusU = ctf.getDefocusU()
            #    defocusV = ctf.getDefocusV()
            #    astigm = defocusU - defocusV
            #    resol = ctf._ctffind4_ctfResolution.get()
            #    if defocusU > self.minDefocus and defocusU < self.maxDefocus and \
            #                defocusV > self.minDefocus and defocusV < self.maxDefocus and \
            #                astigm < self.astigmatism and resol < self.resolution:
            #        outputCTFs.append(ctf)
            #
            #self._defineOutputs(outputCTF=outputCTFs)
            #self._defineTransformRelation(self.ctfSet.get(), outputCTFs)

    def _stepsCheck(self):
        print "en _stepsCheck"
        # Check if there are new micrographs to process
        #ctfSet = self.inputCTFs.get()

        ctfFn = self.inputCTFs.get().getFileName()
        ctfSet = SetOfCTF(filename=ctfFn)
        ctfSet.copyInfo(self.inputCTFs)
        #ctfSet.setMicrographs(self.inputCTFs.get().getMicrographs())

        # Input ctf set can be loaded or None when checked for new inputs
        # If None, we load it
        outputStep = self._getFirstJoinStep()
        newCTFs = self._checkNewCTFs(ctfSet, outputStep)

        if newCTFs is None:
            return

        if outputStep and outputStep.isWaiting():
            outputStep.setStatus(STATUS_NEW)

        #ctfSet.close()


    #--------------------------- INFO functions --------------------------------
    def _summary(self):
        message = []
        return message
    
    def _methods(self):
        pass#nothing here
    
    def _validate(self):
        """ The function of this hook is to add some validation before the protocol
        is launched to be executed. It should return a list of errors. If the list is
        empty the protocol can be executed.
        """
        message = [ ]
        fileCTF = self.inputCTFs.get()
        if fileCTF == None:
            message.append("You must specify a set of CTFs.")
        return message

