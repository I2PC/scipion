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


import pyworkflow.protocol.params as params
import pyworkflow.em as em
from pyworkflow.utils.path import removeBaseExt


class ProtCTFAssign(em.ProtCTFMicrographs):
    """ This protocol assign a CTF estimation to a particular
    set of particles producing a new set. """
    _label = 'ctf assign'
    
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', params.PointerParam, pointerClass='SetOfParticles',
                      label='Input particles',
                      help='Select the particles that you want to update the CTF parameters.')        
        form.addParam('inputCTF', params.PointerParam, pointerClass='SetOfCTF',
                      label="Input CTF",
                      help='Select the CTF that will be used to update particles.')  
        form.addParam('useMicName', params.BooleanParam, default=False,
                      help='Instead of using micId, use the particle stack name\n'
                           'as the link to the correct micrograph. This option\n'
                           'can be useful if the id\'s have been lost due to union.')      

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
        
        defaultCTF = inputCTF.getFirstItem().clone()
        defaultCTF.setDefocusAngle(-9999.)
        defaultCTF.setDefocusU(-9999.)
        defaultCTF.setDefocusV(-9999.)
        
        ctfDict = {}
        for ctf in inputCTF:
            ctfId = ctf.getObjId()
            ctfName = removeBaseExt(ctf.getMicrograph().getFileName())
            ctfDict[ctfId] = ctf
            ctfDict[ctfName] = ctf
            
        missingSet = set() # Report missing micrographs only once
        
        for particle in inputParticles:
            if self.useMicName:
                micKey = particle.getMicId()
            else:
                micKey = removeBaseExt(particle.getFileName())
                
            if micKey not in missingSet:                
                ctf = ctfDict.get(micKey, None)
                
                if ctf is None:
                    self.warning("Discarding particles from micrograph: %s, CTF not found. " % micKey)
                    missingSet.add(micKey)
                else:
                    newParticle = particle.clone()
                    newParticle.setCTF(ctf)
                    outputParticles.append(newParticle)
        
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
