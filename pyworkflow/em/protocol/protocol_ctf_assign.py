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
from pyworkflow.em.protocol import ProtCTFMicrographs


class ProtCTFAssign(ProtCTFMicrographs):
    """ This protocol assign a CTF estimation to a particular
    set of particles producing a new set. """
    _label = 'ctf assign'
    
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles',
                      label='Input particles',
                      help='Select the particles that you want to update the CTF parameters.')        
        form.addParam('inputCTF', PointerParam, pointerClass='SetOfCTF',
                      label="Input CTF",
                      help='Select the CTF that will be used to update particles.')        

        form.addParallelSection(threads=0, mpi=0) 
              
#--------------------------- INSERT steps functions --------------------------------------------  
                                
    def _insertAllSteps(self):
        """for each ctf insert the steps to compare it
        """
        self._insertFunctionStep('createOutputStep')

    def createOutputStep(self):
        inputParticles = self.inputParticles.get()
        inputCTF = self.inputCTF.get() 
        outputParticles = self._createSetOfParticles()
        outputParticles.copyInfo(inputParticles)
        
        for particle in inputParticles:
            micId = particle.getMicId()
            ctf = inputCTF[micId]
            particle.setCTF(ctf)
            outputParticles.append(particle)
        
        self._defineOutputs(outputParticles=outputParticles)
        self._defineSourceRelation(inputParticles, outputParticles)
        self._defineSourceRelation(inputCTF, outputParticles)
        
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
