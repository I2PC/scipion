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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************


from pyworkflow.protocol.params import PointerParam, FloatParam, LEVEL_ADVANCED
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.data import Coordinate
from collections import OrderedDict
from pyworkflow.em.protocol.protocol_particles import ProtParticlePicking

class ProtExtractCoords(ProtParticlePicking):
    """ 
    Extract the coordinates information from a set of particles.
    
    This protocol is useful when we want to re-extract the particles
    (maybe resulting from classification and cleannign) with the 
    original dimensions. It can be also handy to visualize the resulting
    particles in their location in the micrographs.
    """
    _label = 'extract coordinates'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles',
                      label='Input particles', important=True,
                      help='Select the particles from which you want\n'
                           'to extract the coordinates and micrographs.')
        
        form.addParam('inputMicrographs', PointerParam, pointerClass='SetOfMicrographs',
                      label='Input micrographs', important=True,
                      help='Select the micrographs to which you want to\n'
                           'associate the coordinates from the particles.')
        
        form.addParam('correction', FloatParam, level=LEVEL_ADVANCED,
                      label='Scale correction',
                      help='This parameter should not be used. '
                           'It is only now to correct for a wrong '
                           'tracking of the coordinates according to'
                           'the particles pixel size.')

        form.addParallelSection(threads=0, mpi=0)

#--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    def createOutputStep(self):
        inputParticles = self.inputParticles.get()
        inputMics = self.inputMicrographs.get()
        outputCoords = self._createSetOfCoordinates(inputMics)

        scale = inputParticles.getSamplingRate() / inputMics.getSamplingRate()
        c =  self.correction.get() 
        print "Scaling pixel size by: ", scale        
        newCoord = Coordinate()
        
        firstCoord = inputParticles.getFirstItem().getCoordinate()
        hasMicName = firstCoord.getMicName() is not None
            
        # Create the micrographs dict using either micName or micId
        micDict = {}
        missedDict = OrderedDict()
        
        for mic in inputMics:
            micKey = mic.getMicName() if hasMicName else mic.getObjId()
            micDict[micKey] = mic.clone()
            
        for particle in inputParticles:
            coord = particle.getCoordinate()
            micKey = coord.getMicName() if hasMicName else coord.getMicId()
            mic = micDict.get(micKey, None)  
            
            if mic is None: 
                print "Skipping particle, key %s not found" % micKey
            else:
                newCoord.copyObjId(particle)
                #FIXME: the 'correction' as a temporarly hack for fixing the 
                # coordinates scale according to the particles pixel size
                newCoord.setX(coord.getX() * scale * c)
                newCoord.setY(coord.getY() * scale * c)
                newCoord.setMicrograph(mic)
                outputCoords.append(newCoord)
        
        boxSize = inputParticles.getXDim() * scale
        outputCoords.setBoxSize(boxSize)
        
        self._defineOutputs(outputCoordinates=outputCoords)
        self._defineSourceRelation(self.inputParticles, outputCoords)
        self._defineSourceRelation(self.inputMicrographs, outputCoords)
                
    def _summary(self):
        summary = []
        ps1 = self.inputParticles.get().getSamplingRate()
        ps2 = self.inputMicrographs.get().getSamplingRate()
        summary.append('Input particles pixel size: *%0.3f* (A/px)' % ps1)
        summary.append('Input micrographs pixel size: *%0.3f* (A/px)' % ps2)
        summary.append('Scaling coordinates by a factor of *%0.3f*' % (ps1/ps2))
        
        if hasattr(self, 'outputCoordinates'):
            summary.append('Output coordinates: *%d*' % self.outputCoordinates.getSize())
            
        return summary 

    def _methods(self):
        # No much to add to summary information
        return self._summary()

    def _validate(self):
        """ The function of this hook is to add some validation before the protocol
        is launched to be executed. It should return a list of errors. If the list is
        empty the protocol can be executed.
        """
        #check that input set of aligned particles do have 2D alignment
        errors = [ ]
        inputParticles = self.inputParticles.get()
        first = inputParticles.getFirstItem()
        if first.getCoordinate() is None:
            errors.append('The input particles does not have coordinates!!!')
        
        return errors
    
