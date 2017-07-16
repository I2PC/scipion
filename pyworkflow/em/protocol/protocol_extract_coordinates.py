# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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


import numpy
from pyworkflow.protocol.params import PointerParam, BooleanParam
from pyworkflow.em.constants import ALIGN_2D, ALIGN_3D, ALIGN_PROJ, ALIGN_NONE
from pyworkflow.em.data import Coordinate
from pyworkflow.em.protocol.protocol_particles import ProtParticlePicking


class ProtExtractCoords(ProtParticlePicking):
    """ 
    Extract the coordinates information from a set of particles.
    
    This protocol is useful when we want to re-extract the particles
    (maybe resulting from classification and cleaning) with the 
    original dimensions. It can be also handy to visualize the resulting
    particles in their location on micrographs.
    """
    # TESTS: 
    # scipion test tests.em.protocols.test_protocols_xmipp_mics.TestXmippExtractParticles
    
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
        self._insertFunctionStep('createOutputStep')

    def createOutputStep(self):
        inputParticles = self.inputParticles.get()
        inputMics = self.inputMicrographs.get()
        outputCoords = self._createSetOfCoordinates(inputMics)
        alignType = inputParticles.getAlignment()

        scale = inputParticles.getSamplingRate() / inputMics.getSamplingRate()
        print "Scaling coordinates by a factor *%0.2f*" % scale
        newCoord = Coordinate()
        firstCoord = inputParticles.getFirstItem().getCoordinate()
        hasMicName = firstCoord.getMicName() is not None

        # Create the micrographs dict using either micName or micId
        micDict = {}

        for mic in inputMics:
            micKey = mic.getMicName() if hasMicName else mic.getObjId()
            if micKey in micDict:
                print ">>> ERROR: micrograph key %s is duplicated!!!" % micKey
                print "           Used in micrographs:"
                print "           - %s" % micDict[micKey].getLocation()
                print "           - %s" % mic.getLocation()
                raise Exception("Micrograph key %s is duplicated!!!" % micKey)
            micDict[micKey] = mic.clone()

        for particle in inputParticles:
            coord = particle.getCoordinate()
            micKey = coord.getMicName() if hasMicName else coord.getMicId()
            mic = micDict.get(micKey, None)  
            
            if mic is None: 
                print "Skipping particle, key %s not found" % micKey
            else:
                newCoord.copyObjId(particle)
                x, y = coord.getPosition()
                if self.applyShifts:
                    shifts = self.getShifts(particle.getTransform(), alignType)
                    xCoor, yCoor = x - int(shifts[0]), y - int(shifts[1])
                    newCoord.setPosition(xCoor*scale, yCoor*scale)
                else:
                    newCoord.setPosition(x*scale, y*scale)

                newCoord.setMicrograph(mic)
                outputCoords.append(newCoord)
        
        boxSize = inputParticles.getXDim() * scale
        outputCoords.setBoxSize(boxSize)
        
        self._defineOutputs(outputCoordinates=outputCoords)
        self._defineSourceRelation(self.inputParticles, outputCoords)
        self._defineSourceRelation(self.inputMicrographs, outputCoords)

    #--------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        ps1 = self.inputParticles.get().getSamplingRate()
        ps2 = self.inputMicrographs.get().getSamplingRate()
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
        errors = [ ]
        inputParticles = self.inputParticles.get()
        first = inputParticles.getFirstItem()
        if first.getCoordinate() is None:
            errors.append('The input particles do not have coordinates!!!')

        if self.applyShifts and not inputParticles.hasAlignment():
            errors.append('Input particles do not have alignment information!')
        
        return errors

    #--------------------------- UTILS functions ------------------------------
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
        from pyworkflow.em.transformations import translation_from_matrix
        if inverseTransform:
            matrix = numpy.linalg.inv(matrix)
            shifts = -translation_from_matrix(matrix)
        else:
            shifts = translation_from_matrix(matrix)
        return shifts
