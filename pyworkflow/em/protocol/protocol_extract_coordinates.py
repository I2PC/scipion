# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              David Maluenda (dmaluenda@cnb.csic.es)  -streaming version-
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
import numpy
from datetime import datetime

from pyworkflow.protocol.params import PointerParam, BooleanParam
from pyworkflow.protocol.constants import STATUS_NEW
from pyworkflow.em.constants import ALIGN_2D, ALIGN_3D, ALIGN_PROJ, ALIGN_NONE
from pyworkflow.em.data import (Coordinate, SetOfParticles, SetOfMicrographs,
                                SetOfCoordinates, Set)
from pyworkflow.em.protocol import ProtParticlePickingAuto
import pyworkflow.utils as pwutils

import time
class ProtExtractCoords(ProtParticlePickingAuto):
    """ 
    Extract the coordinates information from a set of particles.
    
    This protocol is useful when we want to re-extract the particles
    (maybe resulting from classification and cleaning) with the 
    original dimensions. It can be also handy to visualize the resulting
    particles in their location on micrographs.
    """
    
    _label = 'extract coordinates'

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label='Input particles', important=True,
                      help='Select the particles from which you want\n'
                           'to extract the coordinates and micrographs.')
        
        form.addParam('inputMicrographs', PointerParam,
                      pointerClass='SetOfMicrographs',
                      label='Input micrographs', important=True,
                      help='Select the micrographs to which you want to\n'
                           'associate the coordinates from the particles.')

        form.addParam('applyShifts', BooleanParam, default=False,
                      label='Apply particle shifts?',
                      help='Apply particle shifts from 2D alignment to '
                           'recalculate new coordinates. This can be useful '
                           'for re-centering particle coordinates.')
        
        form.addParallelSection(threads=0, mpi=0)

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self.streamingModeOn = self.getInputParticles().isStreamOpen()

        if self.streamingModeOn:
            self.inputSize = 0
            self.outputSize = 0
            self.micsDone = []

            t0 = time.time()
            newParts, self.streamClosed = self.loadInputs()
            print("loadInputs() time: %fs" % (time.time()-t0))

            stepsIds = self._insertNewSteps(newParts)
            self._insertFunctionStep('createOutputStep',
                                     prerequisites=stepsIds, wait=True)
        else:
            self._insertFunctionStep('createOutputStep')

    def _insertNewSteps(self, partIds):
        deps = []
        stepId = self._insertFunctionStep('extractCoordsStep', partIds,
                                          prerequisites=[])
        deps.append(stepId)
        return deps

    def _stepsCheck(self):
        # To allow streaming picking we need to detect:
        #   1) new micrographs ready to be picked
        #   2) new output coordinates that have been produced and add then
        #      to the output set.
        if self.streamingModeOn:
            self._checkNewInput()
            self._checkNewOutput()
        else:
            pass

    def _checkNewInput(self):
        # Check if there are new particles to process from the input set
        partsFile = self.getInputParticles().getFileName()
        micsFile = self.getInputMicrographs().getFileName()
        now = datetime.now()
        self.lastCheck = getattr(self, 'lastCheck', now)
        mTimeParts = datetime.fromtimestamp(os.path.getmtime(partsFile))
        mTimeMics = datetime.fromtimestamp(os.path.getmtime(micsFile))
        # If the input movies.sqlite have not changed since our last check,
        # it does not make sense to check for new input data
        if self.lastCheck > mTimeParts and self.lastCheck > mTimeMics:
            return None
        self.lastCheck = now

        newParts, self.streamClosed = self.loadInputs()

        if len(newParts) > 0:
            fDeps = self._insertNewSteps(newParts)
            outputStep = self._getFirstJoinStep()
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()

    def extractCoordsStep(self, partsIds):
        outputCoords = self.extractCoordinates(partsIds)

        self.outputSize += len(outputCoords)
        t0 = time.time()
        outputCoords.write()
        outputCoords.close()
        print("write time: %fs" % (time.time()-t0))


    def extractCoordinates(self, partsIds=None):
        inPart = self.getInputParticles()
        inMics = self.getInputMicrographs()
        scale = inPart.getSamplingRate() / inMics.getSamplingRate()
        print "Scaling coordinates by a factor *%0.2f*" % scale
        alignType = inPart.getAlignment()

        suffix = self.getSuffix(partsIds[0]) if partsIds is not None else ''
        outputCoords = self._createSetOfCoordinates(inMics, suffix=suffix)

        def appendCoordFromParticle(part):
            coord = part.getCoordinate()
            micKey = coord.getMicId()
            mic = inMics[micKey]

            if mic is None:
                print "Skipping particle, key %s not found" % micKey
            else:
                newCoord.copyObjId(part)
                x, y = coord.getPosition()
                if self.applyShifts:
                    shifts = self.getShifts(part.getTransform(), alignType)
                    xCoor, yCoor = x - int(shifts[0]), y - int(shifts[1])
                    newCoord.setPosition(xCoor * scale, yCoor * scale)
                else:
                    newCoord.setPosition(x * scale, y * scale)

                newCoord.setMicrograph(mic)
                outputCoords.append(newCoord)

        newCoord = Coordinate()
        if self.streamingModeOn:
            for partId in partsIds:
                particle = inPart[partId]
                appendCoordFromParticle(particle)
        else:
            for particle in inPart:
                appendCoordFromParticle(particle)


        boxSize = inPart.getXDim() * scale
        outputCoords.setBoxSize(boxSize)

        return outputCoords


    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return

        self.finished = self.streamClosed and \
                        self.inputSize == self.outputSize

        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        # we will read all ready files
        files = pwutils.glob(self.getTmpOutputPath('*'))
        newData = len(files) > 0
        lastToClose = self.finished and hasattr(self, 'outputCoordinates')
        if newData or lastToClose:
            outSet = self._loadOutputSet()
            if newData:
                for tmpFile in files:
                    tmpSet = SetOfCoordinates(filename=tmpFile)
                    tmpSet.loadAllProperties()
                    outSet.copyItems(tmpSet)
                    outSet.setBoxSize(tmpSet.getBoxSize())

                    tmpSet.close()
                    pwutils.cleanPath(tmpFile)

            self._updateOutputSet('outputCoordinates', outSet, state=streamMode)

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)

    def _loadOutputSet(self):
        setFile = self._getPath("coordinates.sqlite")
        if os.path.exists(setFile):
            outputSet = SetOfCoordinates(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetOfCoordinates(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)
            self._store(outputSet)
            self._defineTransformRelation(self.getInputParticles(), outputSet)
            self._defineSourceRelation(self.getInputMicrographs(), outputSet)

        outputSet.copyInfo(self.getInputParticles())
        outputSet.setMicrographs(self.getInputMicrographs())

        return outputSet

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None

    def _getFirstJoinStepName(self):
        # This function will be used for streaming, to check which is
        # the first function that need to wait for all micrographs
        # to have completed, this can be overriden in subclasses
        # (e.g., in Xmipp 'sortPSDStep')
        return 'createOutputStep'

    def createOutputStep(self):
        if not self.streamingModeOn:
            outputCoords = self.extractCoordinates()
            self._defineOutputs(outputCoordinates=outputCoords)
            self._defineSourceRelation(self.inputParticles, outputCoords)
            self._defineSourceRelation(self.inputMicrographs, outputCoords)

    # ------------- UTILS functions ----------------
    def getSuffix(self, suffix):
        return "_tmp%s" % suffix

    def getTmpOutputPath(self, suffix):
        return self._getPath("coordinates%s.sqlite" % self.getSuffix(suffix))

    def loadInputs(self):
        micsFn = self.getInputMicrographs().getFileName()
        micsSet = SetOfMicrographs(filename=micsFn)
        micsSet.loadAllProperties()

        availableMics = []
        for mic in micsSet:
            availableMics.append(mic.getObjId())

        micsSetClosed = micsSet.isStreamClosed()
        micsSet.close()

        partsFn = self.getInputParticles().getFileName()
        partsSet = SetOfParticles(filename=partsFn)
        partsSet.loadAllProperties()

        newParts = []
        newMics = []
        for item in partsSet:
            micKey = item.getCoordinate().getMicId()
            if micKey not in self.micsDone and micKey in availableMics:
                newParts.append(item.getObjId())
                if not micKey in self.micsDone:
                    newMics.append(micKey)
        self.micsDone += newMics
        self.inputSize = partsSet.getSize()
        partSetClosed = partsSet.isStreamClosed()
        partsSet.close()

        return newParts, micsSetClosed and partSetClosed

    def getShifts(self, transform, alignType):
        """
        is2D == True-> matrix is 2D (2D images alignment)
                otherwise matrix is 3D (3D volume alignment or projection)
        invTransform == True  -> for xmipp implies projection
                              -> for xmipp implies alignment
        """
        if alignType == ALIGN_NONE:
            return None

        inverseTransform = alignType == ALIGN_PROJ
        # only flip is meaningful if 2D case
        # in that case the 2x2 determinant is negative
        flip = False
        matrix = transform.getMatrix()
        if alignType == ALIGN_2D:
            # get 2x2 matrix and check if negative
            flip = bool(numpy.linalg.det(matrix[0:2, 0:2]) < 0)
            if flip:
                matrix[0, :2] *= -1.  # invert only the first two columns keep x
                matrix[2, 2] = 1.  # set 3D rot
            else:
                pass

        elif alignType == ALIGN_3D:
            flip = bool(numpy.linalg.det(matrix[0:3, 0:3]) < 0)
            if flip:
                matrix[0, :4] *= -1.  # now, invert first line including x
                matrix[3, 3] = 1.  # set 3D rot
            else:
                pass

        else:
            pass
            # flip = bool(numpy.linalg.det(matrix[0:3,0:3]) < 0)
            # if flip:
            #    matrix[0,:4] *= -1.#now, invert first line including x
        shifts = self.geometryFromMatrix(matrix, inverseTransform)

        return shifts

    def geometryFromMatrix(self, matrix, inverseTransform):
        from pyworkflow.em.convert.transformations import translation_from_matrix
        if inverseTransform:
            matrix = numpy.linalg.inv(matrix)
            shifts = -translation_from_matrix(matrix)
        else:
            shifts = translation_from_matrix(matrix)
        return shifts

    def getInputParticles(self):
        return self.inputParticles.get()

    #--------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        ps1 = self.getInputParticles().getSamplingRate()
        ps2 = self.getInputMicrographs().getSamplingRate()
        summary.append(u'Input particles pixel size: *%0.3f* (Å/px)' % ps1)
        summary.append(u'Input micrographs pixel size: *%0.3f* (Å/px)' % ps2)
        summary.append('Scaling coordinates by a factor of *%0.3f*' % (ps1/ps2))
        if self.applyShifts:
            summary.append('Applied 2D shifts from particles')

        if hasattr(self, 'outputCoordinates'):
            summary.append('Output coordinates: *%d*'
                           % self.outputCoordinates.getSize())

        return summary

    def _methods(self):
        # No much to add to summary information
        return self._summary()

    def _validate(self):
        """ The function of this hook is to add some validation before the
        protocol is launched to be executed. It should return a list of errors.
        If the list is empty the protocol can be executed.
        """
        errors = []
        inputParticles = self.getInputParticles()
        first = inputParticles.getFirstItem()
        if first.getCoordinate() is None:
            errors.append('The input particles do not have coordinates!!!')

        if self.applyShifts and not inputParticles.hasAlignment():
            errors.append('Input particles do not have alignment information!')

        return errors
