# **************************************************************************
# *
# * Authors:     Laura del Cano (laura.cano@cnb.csic.es)
# *              Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *              I. Foche (ifoche@cnb.csic.es)
# *              Tomas Majtner (tmajtner@cnb.csic.es)   -- streaming version
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
import pyworkflow.em.metadata as md
import pyworkflow.protocol.constants as cons
from pyworkflow.em.data import SetOfParticles
from pyworkflow.em.protocol import ProtProcessParticles
from pyworkflow.object import String, Set
from pyworkflow.protocol.params import (EnumParam, IntParam, Positive, Range,
                                        LEVEL_ADVANCED, FloatParam, BooleanParam)
from convert import readSetOfParticles, writeSetOfParticles, setXmippAttributes

class XmippProtScreenParticles(ProtProcessParticles):
    """ Classify particles according their similarity to the others in order
    to detect outliers. """

    _label = 'screen particles'

    # Automatic Particle rejection enum
    REJ_NONE = 0
    REJ_MAXZSCORE = 1
    REJ_PERCENTAGE =2
    REJ_PERCENTAGE_SSNR =1
    #--------------------------- DEFINE param functions ----------------------
    def _defineProcessParams(self, form):
        
        form.addParam('autoParRejection', EnumParam,
                      choices=['None', 'MaxZscore', 'Percentage'],
                      label="Automatic particle rejection based on Zscore",
                      default=self.REJ_NONE,
                      display=EnumParam.DISPLAY_COMBO,
                      expertLevel=LEVEL_ADVANCED,
                      help='How to automatically reject particles. It can be:\n'
                           '  None (no rejection)\n'
                           '  MaxZscore (reject a particle if its Zscore [a '
                           'similarity index] is larger than this value).\n '
                           '  Percentage (reject a given percentage in each '
                           'one of the screening criteria).')
        form.addParam('maxZscore', FloatParam, default=3,
                      condition='autoParRejection==1',
                      label='Maximum Zscore', expertLevel=LEVEL_ADVANCED,
                      help='Maximum Zscore.', validators=[Positive])      
        form.addParam('percentage', IntParam, default=5,
                      condition='autoParRejection==2',
                      label='Percentage (%)', expertLevel=LEVEL_ADVANCED,
                      help='The worse percentage of particles according to '
                           'metadata labels: ZScoreShape1, ZScoreShape2, '
                           'ZScoreSNR1, ZScoreSNR2, ZScoreHistogram are '
                           'automatically disabled. Therefore, the total '
                           'number of disabled particles belongs to ['
                           'percetage, 5*percentage]',
                      validators=[Range(0, 100, error="Percentage must be "
                                                      "between 0 and 100.")])
        form.addParam('autoParRejectionSSNR', EnumParam,
                      choices=['None', 'Percentage'],
                      label="Automatic particle rejection based on SSNR",
                      default=self.REJ_NONE, display=EnumParam.DISPLAY_COMBO,
                      expertLevel=LEVEL_ADVANCED,
                      help='How to automatically reject particles based on '
                           'SSNR. It can be:\n'
                           '  None (no rejection)\n'
                           'Percentage (reject a given percentage of the '
                           'lowest SSNRs).')
        form.addParam('percentageSSNR', IntParam, default=5,
                      condition='autoParRejectionSSNR==1',
                      label='Percentage (%)', expertLevel=LEVEL_ADVANCED,
                      help='The worse percentage of particles according to '
                           'SSNR are automatically disabled.',
                      validators=[Range(0, 100, error="Percentage must be "
                                                      "between 0 and 100.")])
        form.addParam('addFeatures', BooleanParam, default=False,
                      label='Add features', expertLevel=LEVEL_ADVANCED,
                      help='Add features used for the ranking to each one '
                           'of the input particles')
        form.addParallelSection(threads=0, mpi=0)
        
    def _getDefaultParallel(self):
        """This protocol doesn't have mpi version"""
        return (0, 0)
     
    #--------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self.insertedDict = {}
        self.finished = False
        self.SetOfParticles = [m.clone() for m in self.inputParticles.get()]
        partsSteps = self._insertNewPartsSteps(self.insertedDict,
                                               self.SetOfParticles)
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
        writeSetOfParticles([p.clone() for p in inputParts
                             if int(p.getObjId()) not in insertedDict],
                            self._getPath('images.xmd'),
                            alignType=em.ALIGN_NONE)
        stepId = self._insertFunctionStep('sortImages',
                                          self._getPath('images.xmd'),
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
        self.finished = self.streamClosed and (len(doneList) == len(partsSet))
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

    def _loadOutputSet(self, SetClass, baseName, fnInputMd):

        setFile = self._getPath(baseName)
        if os.path.exists(setFile):
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        partsSet = self._createSetOfParticles()
        readSetOfParticles(fnInputMd, partsSet)
        inputs = self.inputParticles.get()
        outputSet.copyInfo(inputs)
        outputSet.copyItems(partsSet,
                            updateItemCallback=self._updateParticle,
                            itemDataIterator=md.iterRows(
                                self.outputMd.get(),
                                sortByLabel=md.MDL_ITEM_ID))
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

    #--------------------------- STEPS functions -----------------------------
    def sortImages(self, fnInputMd):

        args = "-i Particles@%s --addToInput " % fnInputMd
        
        if self.autoParRejection == self.REJ_MAXZSCORE:
            args += "--zcut " + str(self.maxZscore.get())
        
        elif self.autoParRejection == self.REJ_PERCENTAGE:
            args += "--percent " + str(self.percentage.get())

        if self.addFeatures:
            args += "--addFeatures "

        self.runJob("xmipp_image_sort_by_statistics", args)

        self.outputMd = String(fnInputMd)

        args = "-i Particles@%s " % fnInputMd
        
        if self.autoParRejectionSSNR == self.REJ_PERCENTAGE_SSNR:
            args += "--ssnrpercent " + str(self.percentageSSNR.get())

        self.runJob("xmipp_image_ssnr", args)

        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN
        outSet = self._loadOutputSet(SetOfParticles, 'outputParticles.sqlite',
                                     fnInputMd)
        self._updateOutputSet('outputParticles', outSet, streamMode)

    #--------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        else:
            fnSummary = self._getExtraPath("summary.txt")
            if not os.path.exists(fnSummary):
                zscores = [p._xmipp_zScore.get() for p in self.outputParticles]
                if len(zscores)>0:
                    fhSummary = open(fnSummary,"w")
                    fhSummary.write("The minimum ZScore is %.2f\n" % min(zscores))
                    fhSummary.write("The maximum ZScore is %.2f\n" % max(zscores))
                    fhSummary.write("The mean ZScore is %.2f\n"
                                    % (sum(zscores)*1.0/len(self.outputParticles)))
                fhSummary.close()
            if os.path.exists(fnSummary):
                fhSummary = open(fnSummary)
                for line in fhSummary.readlines():
                    summary.append(line.strip())
                fhSummary.close()
        return summary
    
    def _validate(self):
        pass
        
    def _citations(self):
        return ['Vargas2013b']
    
    def _methods(self):
        methods = []
        if hasattr(self, 'outputParticles'):
            outParticles = (len(self.outputParticles) if self.outputParticles
                                                         is not None else None)
            particlesRejected = (len(self.inputParticles.get())-outParticles
                                 if outParticles is not None else None)
            particlesRejectedText = (' ('+str(particlesRejected)+')' if
                                     particlesRejected is not None else '')
            rejectionText = ['',# REJ_NONE
                             ' and removing those not reaching %s%s'
                             % (str(self.maxZscore.get()),
                                particlesRejectedText),# REJ_MAXZSCORE
                             ' and removing worst %s percent%s'
                             % (str(self.percentage.get()),
                                particlesRejectedText)# REJ_PERCENTAGE
                             ]
            methods.append('Input dataset %s of %s particles was sorted by'
                           ' its ZScore using xmipp_image_sort_by_statistics'
                           ' program%s. '
                           % (self.getObjectTag('inputParticles'),
                              len(self.inputParticles.get()),
                              rejectionText[self.autoParRejection.get()]))
            methods.append('Output set is %s.'
                           % self.getObjectTag('outputParticles'))
        return methods
    
    #--------------------------- UTILS functions -----------------------------
    def _updateParticle(self, item, row):
        setXmippAttributes(item, row, md.MDL_ZSCORE, md.MDL_ZSCORE_SHAPE1,
                           md.MDL_ZSCORE_SHAPE2, md.MDL_ZSCORE_SNR1,
                           md.MDL_ZSCORE_SNR2, md.MDL_CUMULATIVE_SSNR)
        if self.addFeatures:
            setXmippAttributes(item, row, md.MDL_SCORE_BY_SCREENING)
        if row.getValue(md.MDL_ENABLED) <= 0:
            item._appendItem = False
        else:
            item._appendItem = True

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