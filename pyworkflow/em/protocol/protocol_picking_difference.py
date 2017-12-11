# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
from pyworkflow.em.protocol.protocol_particles import ProtParticlePicking


class ProtPickingDifference(ProtParticlePicking):
    """
    Protocol to compute the difference between a reference SetOfPartices and
    a another set (usually a negative reference).

    The output will be a SetOfCoordinates with the particles in the reference
    input that are not close to coordinates in the negative input.

    """
    _label = 'picking difference'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates',
                      label="Input coordinates",
                      help='Select the reference set of coordinates.')
        form.addParam('negativeCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates',
                      label="Negative coordinates",
                      help='Negative coordinates that will help to exclude '
                           'coordinates from the input set. ')
        form.addParam('differenceRadius',params.IntParam, default=10,
                      label="Radius (px)",
                      help="Distance radius to consider that two coordinates "
                           "close enough. If a particle in the input reference"
                           "set have any close particle in the negative set, "
                           "it will not be included in the output set. ")

        form.addParallelSection(threads=0, mpi=0)
        
    # ------------------------- INSERT steps functions -------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep("createOutputStep",
                                 self.inputCoordinates.getObjId(),
                                 self.negativeCoordinates.getObjId(),
                                 self.differenceRadius.get())

    def getInputMicrographs(self):
        return self.inputCoordinates[0].get().getMicrographs()
    
    def _summary(self):
        summary = []
        return summary
    
    def _methods(self):
        return []
    
    def createOutputStep(self, inputId, negId, radius):
        inputCoords = self.inputCoordinates.get()
        negCoords = self.negativeCoordinates.get()
        inputMics = inputCoords.getMicrographs()
        outputCoords = self._createSetOfCoordinates(inputMics)
        outputCoords.setBoxSize(inputCoords.getBoxSize())

        r2 = radius * radius # to avoid computing sqrt when comparing distances

        def _discard(coord, negPos):
            """ Find if there is a close negative particle to this one. """
            x, y = coord.getPosition()

            for nx, ny in negPos:
                if abs((x - nx) * (y - ny)) < r2:
                    return True

            return False # Far enough from all negative coordinates

        # Read all consensus particles
        for mic in inputMics:
            negPos = [(c.getX(), c.getY())
                      for c in negCoords.iterCoordinates(mic)]
            for coord in inputCoords.iterCoordinates(mic):
                if not _discard(coord, negPos):
                    outputCoords.append(coord)

        # Set output
        self._defineOutputs(outputCoordinates=outputCoords)
        self._defineTransformRelation(self.inputCoordinates,
                                      self.outputCoordinates)
        self._defineSourceRelation(self.negativeCoordinates,
                                   self.outputCoordinates)
