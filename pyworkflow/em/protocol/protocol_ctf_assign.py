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
from pyworkflow.em.data import SetOfParticles
from pyworkflow.em.protocol import ProtCTFMicrographs


class ProtCTFAssign(ProtCTFMicrographs):
    """ This protocol assign a CTF estimation to a particular
    set of particles producing a new set. """
    _label = 'ctf assign'
    _unionTypes = ['Micrographs',
                   'Particles']
    
    def __init__(self, **kwargs):
        ProtCTFMicrographs.__init__(self, **kwargs)
        # We need to trace the changes of 'inputType' to 
        # dynamically modify the property of pointerClass
        # of the 'inputSets' parameter
        def onChangeInputType():
            pointerClass = 'SetOf' + self.getEnumText('inputType')
            self.inputSetsParam.setPointerClass(pointerClass)
        # Initial update
        onChangeInputType()
        # Now keep track of changes and update
        self.inputType.trace(onChangeInputType)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputType', params.EnumParam, choices=self._unionTypes, default=0, # Micrographs
                      label='Input type:',
                      help='Select the type of objects that you want to assign the CTF.')
        self.inputSetsParam = form.addParam('inputSet', params.PointerParam, pointerClass='EMSet',
                              label='Input set',
                              help='Select the images (micrographs or particles) '
                                   'that you want to update the CTF parameters.')        
        form.addParam('inputCTF', params.PointerParam, pointerClass='SetOfCTF',
                      label="Input CTF",
                      help='Select the CTF that will be used to update particles.')  

        form.addParallelSection(threads=0, mpi=0) 
        
#--------------------------- INSERT steps functions --------------------------------------------  
                                
    def _insertAllSteps(self):
        """for each ctf insert the steps to compare it
        """
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self):
        inputSet = self.inputSet.get()
        inputCTF = self.inputCTF.get()
        
        if isinstance(inputSet, SetOfParticles):
            self._particlesOutputStep(inputSet, inputCTF)
        else:
            self._microgrpahsOutputStep(inputSet, inputCTF)
            
    def _particlesOutputStep(self, inputSet, inputCTF):
        outputParts = self._createSetOfParticles()
        outputParts.copyInfo(inputSet)
        
        ctfDict = {}
        for ctf in inputCTF:
            ctfName = ctf.getMicrograph().getMicName()
            ctfDict[ctfName] = ctf
        
        missingSet = set() # Report missing micrographs only once
        
        firstCoord = inputSet.getFirstItem().getCoordinate()
        hasMicName = firstCoord.getMicName() is not None
        
        for particle in inputSet:
            coord = particle.getCoordinate()
            micKey = coord.getMicName() if hasMicName else coord.getObjId()
            
            if micKey not in missingSet:
                ctf = ctfDict.get(micKey, None)
                
                if ctf is None:
                    self.warning("Discarding particles from micrograph with micName: %s, CTF not found. " % micKey)
                    missingSet.add(micKey)
                else:
                    newParticle = particle.clone()
                    newParticle.setCTF(ctf)
                    outputParts.append(newParticle)
        
        self._defineOutputs(outputParticles=outputParts)
        self._defineSourceRelation(self.inputSet, outputParts)
        self._defineSourceRelation(self.inputCTF, outputParts)
    
    def _microgrpahsOutputStep(self, inputSet, inputCTF):           
        outputMics = self._createSetOfMicrographs()
        outputMics.copyInfo(inputSet)
        ctfDict = {}
        
        for ctf in inputCTF:
            ctfName = ctf.getMicrograph().getMicName()
            ctfDict[ctfName] = ctf
        
        for mic in inputSet:
            micKey = mic.getMicName()
            ctf = ctfDict.get(micKey, None)
            if ctf is None:
                self.warning("Discarding micrographs with micName: %s, CTF not found. " % micKey)
            else:
#                 ctf.setObjId(mic.getObjId())
                newMic = mic.clone()
                outputMics.append(newMic)
        
        self._defineOutputs(outputMicrographs=outputMics)
        self._defineSourceRelation(self.inputSet, outputMics)
        self._defineCtfRelation(outputMics, self.inputCTF)
    
    #--------------------------- INFO functions ----------------------------------------------------
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
