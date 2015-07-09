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


from pyworkflow.protocol.params import PointerParam, FloatParam
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.data import Coordinate


class ProtExtractCoords(EMProtocol):
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
        
        form.addParam('correction', FloatParam, 
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
        
        for particle in inputParticles:
            coord = particle.getCoordinate()
            micId = particle.getMicId()
            mic = inputMics[micId]  
            
            if mic is None:
                print "Skipping particle, micId %s not found" % micId
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
        self._defineSourceRelation(inputParticles, outputCoords)
        self._defineSourceRelation(inputMics, outputCoords)
                
    def _summary(self):
        summary = []
        return summary 
    
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        else:
            scale = self.inputAlignment.get().getSamplingRate()/self.inputParticles.get().getSamplingRate()
            summary.append("Assigned alignment to %s particles from a total of %s." % (self.outputParticles.getSize(), self.inputParticles.get().getSize()))
            if scale != 1:
                summary.append("Applied scale of %s." % scale)
        return summary

    def _methods(self):
        methods = []
        return methods
    
        if not hasattr(self, 'outputParticles'):
            methods.append("Output particles not ready yet.")
        else:
            scale = self.inputAlignment.get().getSamplingRate()/self.inputParticles.get().getSamplingRate()
            methods.append("We assigned alignment to %s particles from %s and produced %s."
                           % (self.outputParticles.getSize(), self.getObjectTag('inputParticles'), self.getObjectTag('outputParticles')))
            if scale != 1:
                methods.append("Applied scale factor of %s." % scale)
        return methods

    def _validate(self):
        """ The function of this hook is to add some validation before the protocol
        is launched to be executed. It should return a list of errors. If the list is
        empty the protocol can be executed.
        """
        #check that input set of aligned particles do have 2D alignment
        errors = [ ]
        return errors
    
        inputAlignmentSet = self.inputAlignment.get()
        if not inputAlignmentSet.hasAlignment():
            errors.append("Input alignment set should contains some kind of alignment (2D, 3D or Projection).")
        else:
            # Just for consistency, check that the particles really contains Transform object
            first = inputAlignmentSet.getFirstItem()
            alignment = first.getTransform()
            if alignment is None:
                errors.append('Inconsistency detected in *Input alignment* !!!')
                errors.append('It has alignment: _%s_, but the alignment is missing!!!' % inputAlignmentSet.getAlignment())
            
        # Add some errors if input is not valid
        return errors
