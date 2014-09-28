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


from pyworkflow.protocol.params import PointerParam
from protocol_2d import ProtAlign2D


class ProtAlignmentAssign(ProtAlign2D):
    """ Assign a the alignment calculated for a set of particles to another set.
    This protocol will take into account the differences of pixel size (A/pix)
    between the two sets and multiply by the right factor the shifts.
    The particles with the alignment can also be a subset of the other images
    """
    _label = 'alignment assign'
    
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles',
                      label='Input particles',
                      help='Select the particles that you want to update the new alignment.')        
        form.addParam('inputAlignment', PointerParam, pointerClass='SetOfParticles',
                      label="Input alignments",
                      help='Select the particles with alignment to be apply to the other particles.')        

        form.addParallelSection(threads=0, mpi=0) 
              
#--------------------------- INSERT steps functions --------------------------------------------  
                                
    def _insertAllSteps(self):
        """for each ctf insert the steps to compare it
        """
        self._insertFunctionStep('createOutputStep')

    def createOutputStep(self):
        inputParticles = self.inputParticles.get()
        inputAlignment = self.inputAlignment.get() 
        outputParticles = self._createSetOfParticles()
        outputParticles.copyInfo(inputParticles)
        
        scale = inputAlignment.getSamplingRate()/inputParticles.getSamplingRate()
        n = inputParticles.getSize()
        block = min(n/10, 1000)       
        
        for i, particle in enumerate(inputParticles):
            alignedParticle = inputAlignment[particle.getObjId()]
            if alignedParticle is not None:
                newParticle = particle.clone()
                alignment = alignedParticle.getAlignment()
                alignment._xmipp_shiftX.multiply(scale)
                alignment._xmipp_shiftY.multiply(scale)
                newParticle.setAlignment(alignment)
                outputParticles.append(newParticle)
            ii = i+1
            if ii % block == 0:
                self.info('Done %d out of %d' % (i+1, n))
        
        self._defineOutputs(outputParticles=outputParticles)
        self._defineSourceRelation(inputParticles, outputParticles)
        self._defineSourceRelation(inputAlignment, outputParticles)
        
    def _summary(self):
        summary = []
        return summary    
    
    def _methods(self):
        return []
    
    def _validate(self):
        """ The function of this hook is to add some validation before the protocol
        is launched to be executed. It should return a list of errors. If the list is
        empty the protocol can be executed.
        """
        #same micrographs in both CTF??
        errors = [ ] 
        # Add some errors if input is not valid
        return errors
