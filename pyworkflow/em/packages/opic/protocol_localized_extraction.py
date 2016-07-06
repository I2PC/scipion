# *****************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# *****************************************************************************
"""
This module contains the protocol for localized reconstruction.
"""

import numpy as np

from pyworkflow.em import ALIGN_PROJ, ImageHandler
from pyworkflow.protocol.params import PointerParam
from pyworkflow.em.protocol import ProtParticles, IntParam



class ProtLocalizedExtraction(ProtParticles):
    """ Extract computed sub-particles from a SetOfParticles.
    """
    
    _label = 'localized extraction'
    
    #--------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      important=True,
                      label="Input particles",  
                      help='Select the input images from the project.')

        form.addParam('inputCoordinates', PointerParam,
                      pointerClass='SetOfCoordinates',
                      important=True,
                      label='Input coordinates')

        form.addParam('boxSize', IntParam,
                      label='Subparticle box size (px)',
                      help='Select the amount of pixels to extract the '
                           'sub-particles.')

        form.addParallelSection(threads=0, mpi=0)

    #--------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        partsId = self.inputParticles.get().getObjId()
        self._insertFunctionStep('createOutputStep',
                                 self._getInputParticles().getObjId(),
                                 self.inputCoordinates.get().getObjId(),
                                 self.boxSize.get())
    
    #--------------------------- STEPS functions ------------------------------
    def createOutputStep(self, particlesId, coordsId, boxSize):
        """ Create the input file in STAR format as expected by Relion.
        Params:
            particlesId: use this parameters just to force redo of convert if 
                the input particles are changed.
        """
        ih = ImageHandler()
        outputStack = self._getPath('particles.mrcs')
        outputImg = ih.createImage()

        inputParticles = self._getInputParticles()
        inputCoords = self.inputCoordinates.get()
        outputSet = self._createSetOfParticles()
        outputSet.copyInfo(inputParticles)

        boxSize = self.boxSize.get()
        b2 = int(round(boxSize / 2))
        center = np.zeros((boxSize, boxSize))

        ih = ImageHandler()

        i = 0
        lastPartId = None

        for coord in inputCoords.iterItems(orderBy=['_subparticle._micId',
                                                    '_micId', 'id']):
            # The original particle id is stored in the sub-particle as micId
            partId = coord._micId.get()

            # Load the particle if it has changed from the last sub-particle
            if partId != lastPartId:
                particle = inputParticles[partId]

                if particle is None:
                    raise Exception('Missing particle with id %s from '
                                    'input particles set' % partId)
                lastPartId = partId
                # Now load the particle image to extract later sub-particles
                img = ih.read(particle)
                data = img.getData()

            xpos = coord.getX()
            ypos = coord.getY()

            # Crop the sub-particle data from the whole particle image
            center[:, :] = data[ypos-b2:ypos+b2, xpos-b2:xpos+b2]
            outputImg.setData(center)
            i += 1
            outputImg.write((i, outputStack))
            subpart = coord._subparticle
            subpart.setLocation((i, outputStack)) # Change path to new stack
            subpart.setObjId(None) # Force to insert as a new item
            outputSet.append(subpart)

        self._defineOutputs(outputParticles=outputSet)
        self._defineSourceRelation(self.inputParticles, outputSet)

    #--------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        inputCoords = self.inputCoordinates.get()
        firstCoord = inputCoords.getFirstItem()

        if not firstCoord.hasAttribute('_subparticle'):
            errors.append('The selected input coordinates does not are the'
                          'output from a localized-subparticles protocol.')

        return errors
    
    def _citations(self):
        return ['Serban2015']

    def _summary(self):
        summary = []
        return summary
    
    def _methods(self):
        return []
    
    #--------------------------- UTILS functions ------------------------------
    def _getInputParticles(self):
        return self.inputParticles.get()
    
